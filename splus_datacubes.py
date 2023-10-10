## Main code by Kadu Barbosa
## Changes made by Amanda Lopes

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import astropy.units as u
import os
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
import astropy.constants as const
from astropy.wcs import WCS
from scipy.interpolate import RectBivariateSpline
from tqdm import tqdm
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from matplotlib.colors import LogNorm
from spatial_el import SPLUSEmLine
warnings.simplefilter('ignore', category=AstropyWarning)
np.seterr(divide='ignore', invalid='ignore')

class SCube():
    """ Produces S-PLUS data cubes directly from the project's database.

    Parameters
    ----------
    obj: str
        Identification of the object
    coords: astropy.coordinates.SkyCoord or list
        Coordinates of the object. If list is given, it assumes units in
        degrees for both coordinates.
    size: astropy.Quantity or flot
        Size of the cube.
    survey: str
        Survey, for example S-PLUS, J-PLUS, J-PAS (future updates)
    field: str
        Field of the object
    conn: splusdata.conn
        Connection with S-PLUS data base using splusdata
    wdir: str
        Working directory
    zpref: str
        Calibration reference. Options are idr3_n4 (latest calibration,
        default) and idr3
    redo: bool
        Produces cube again even if it already exists in the disk.
        Default is False.
    bands: list
        Names of the filters to be included in the cube. Options include
        'U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I',
        'F861', and 'Z'. Defalt is all bands.
    coord_unit: astropy.unit
        Units to be used in coordinates. Default is degree
    cut_dir: str
        Path to location where cutouts are stored.
    """
    def __init__(self, obj, coords, size, survey, field=None, conn=None, wdir=None, zpref=None,
                 bands=None, verbose=True, coord_unit=None,
                 cut_dir=None):

        self.obj = obj
        self.coord_unit = u.degree if coord_unit is None else coord_unit
        if isinstance(coords, SkyCoord):
            self.coords = coords
        elif isinstance(coords, list):
                self.coords = SkyCoord(coords[0], coords[1],
                                       unit=self.coord_unit)
        else:
            raise ValueError("Input coordinates should be list or SkyCoord")
        self.size = size
        if isinstance(self.size, u.Quantity):
            self.size_unit = self.size.unit
        else:
            self.size_unit = u.pix
        self.survey = survey
        self.conn = conn
        self.field = field
        if self.field==None:
            self.field = 1
        self.verbose = verbose
        # General definitions
        self.ps = 0.55 * u.arcsec / u.pixel
        if isinstance(self.size, u.Quantity):
            self.cutsize = int((self.size / self.ps).value)
        else:
            self.cutsize = int(self.size)
        if (self.survey == 'S-PLUS') or (self.survey == 'J-PLUS'):
            all_bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R',
                      'F660', 'I', 'F861', 'Z']
            if bands is not None:
                    assert all([x in all_bands for x in bands]), \
                   "Input bands list should include only S-PLUS/J-PLUS filters"
            self.bands_names = {'U': "$u$", 'F378': "$J378$", 'F395': "$J395$",
                            'F410': "$J410$", 'F430': "$J430$", 'G': "$g$",
                            'F515': "$J515$", 'R': "$r$", 'F660': "$J660$",
                            'I': "$i$", 'F861': "$J861$", 'Z': "$z$"}
            self.wave_eff = {"F378": 3770.0, "F395": 3940.0, "F410": 4094.0,
                         "F430": 4292.0, "F515": 5133.0, "F660": 6614.0,
                         "F861": 8611.0, "G": 4751.0, "I": 7690.0, "R": 6258.0,
                         "U": 3536.0, "Z": 8831.0}
        self.bands = all_bands if bands is None else bands
        self.wave = np.array([self.wave_eff[band] for band in
                              self.bands]) * u.Angstrom
        self.flam_unit = u.erg / u.cm / u.cm / u.s / u.AA
        self.fnu_unit = u.erg / u.s / u.cm / u.cm / u.Hz
        # Setting directory for stamps
        self.wdir = os.getcwd() if wdir is None else wdir
        self.cutouts_dir = os.path.join(self.wdir, "cutouts") \
                            if cut_dir is None else cut_dir
        # Setting zero point calibration configurations
        self._path = os.path.dirname(os.path.abspath(__file__))
        self.zpref = "idr3_n4" if zpref is None else zpref
        self.zps = self.get_zps(self.wdir)
        self.zpcorr = self.get_zp_correction(self.wdir)

        if not os.path.exists(self.cutouts_dir):
            os.mkdir(self.cutouts_dir)
        # Producing stamps and cube
        self._s = int(self.size)
        self.cutnames = [f"{self.obj}_{band}_{self._s}x{self._s}" \
                         f"{self.size_unit}.fz" for band in self.bands]
        self.wcutnames= [cut.replace(".fz", "_weight.fz") for cut in
                         self.cutnames]
        self.output_dir = f"{self.wdir}/output/"
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.cubename = os.path.join(self.output_dir, f"{self.obj}_{self._s}x{self._s}"
                                                f"{self.size_unit}.fits")


    def get_status(self):
        """ Check if proposed objects are found in the S-PLUS DB. """
        if self.conn is None:
            return 0, "Connection with SPLUS DB not provided."
        ra = self.coords.ra.to(u.degree).value
        dec = self.coords.dec.to(u.degree).value
        try:
            self.conn.get_cut(ra, dec, 5, self.bands[0])
            return 1, None
        except:
            return 0, "Region not available in S-PLUS database."

    def download_stamps(self, redo=False):
        status, msg = self.get_status()
        if not status:
            print(msg)
            return
        desc = "Downloading stamps"
        N = len(self.bands)
        for band, img, wimg in tqdm(zip(self.bands, self.cutnames,
                                      self.wcutnames), desc=desc, total=N):
            ra = self.coords.ra.to(u.degree).value
            dec = self.coords.dec.to(u.degree).value
            outfile = os.path.join(self.cutouts_dir, img)
            if not os.path.exists(outfile) or redo:
                hdu = self.conn.get_cut(ra, dec, self.cutsize, band, option=self.field)
                hdu.writeto(outfile, overwrite=True)
            # Getting the weight images.
            woutfile = os.path.join(self.cutouts_dir, wimg)
            if not os.path.exists(woutfile) or redo:
                hdu = self.conn.get_cut_weight(ra, dec, self.cutsize, band, option=self.field)
                hdu.writeto(woutfile, overwrite=True)

    def get_zps(self,dirwork):
        """ Load all tables with zero points for iDR3. """
        zp_dir = os.path.join(dirwork, f"assets/{self.zpref}/zps")
        tables = []
        for fname in os.listdir(zp_dir):
            filename = os.path.join(zp_dir, fname)
            data = np.genfromtxt(filename, dtype=None, encoding=None)
            with open(filename) as f:
                h = f.readline().replace("#", "").replace("SPLUS_", "").split()
            table = Table(data, names=h)
            tables.append(table)
        zptable = vstack(tables)
        return zptable

    def get_zp_correction(self,dirwork):
        """ Get corrections of zero points for location in the field. """
        x0, x1, nbins = 0, 9200, 16
        xgrid = np.linspace(x0, x1, nbins + 1)
        zpcorr = {}
        for band in self.bands:
            corrfile = os.path.join(dirwork,
                        f"assets/zpcorr_idr3/SPLUS_{band}_offsets_grid.npy")
            corr = np.load(corrfile)
            zpcorr[band] = RectBivariateSpline(xgrid, xgrid, corr)
        return zpcorr

    def remove_stars(self,data):
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)
        sources = daofind(data - median)
        xdim, ydim = data.shape
        x0 = xdim / 2.
        y0 = ydim / 2.
        r = np.sqrt((sources["xcentroid"] - x0) ** 2 + (sources["ycentroid"] - y0) ** 2).data
        idx = np.where(r >= 90)
        stars = sources[idx]
        x = np.arange(xdim)
        y = np.arange(ydim)
        xx, yy = np.meshgrid(x, y)
        mask_star = np.zeros_like(data).astype(np.bool)
        rstars = 6
        for star in stars:
            r = np.sqrt((xx - star["xcentroid"]) ** 2 + (yy - star["ycentroid"]) ** 2)
            idx = np.where(r < rstars)
            mask_star[idx] = True
        return mask_star

    def make_cube_correted(self, el=None, redo=False):
        dir_corrected = os.path.join(self.cutouts_dir, f"pre_processing/psf_corrected/{self.obj}/")
        if os.path.exists(self.cubename) and not redo:
            return
        cuts_exist = [os.path.exists(os.path.join(self.cutouts_dir, img)) for
                      img in self.cutnames]
        if not all(cuts_exist):
            print("Not all stamps are available locally, cube not produced.")
            return
        fields = [_.replace("_", "-") for _ in self.zps["FIELD"]]
        flams, flamerrs = [], []
        headers = []
        for band, img in zip(self.bands, self.cutnames):
            wave = self.wave_eff[band] * u.Angstrom
            filename = os.path.join(self.cutouts_dir, img)
            h = fits.getheader(filename, ext=1)
            tile = h["OBJECT"].replace("_", "-")
            # Getting zero-point calibration and location correction
            zp0 = self.zps[fields.index(tile)][band]
            x0 = h["X0TILE"]
            y0 = h["Y0TILE"]
            zpcorr = float(self.zpcorr[band](x0, y0))
            zp = zp0 + zpcorr
            gain = h["GAIN"]
            f0 = np.power(10, -0.4 * (48.6 + zp))
            file_projected = os.path.join(dir_corrected, f"{self.obj}_{band}_{self._s}x{self._s}{self.size_unit}.fits")
            data = fits.getdata(file_projected, 0)
            weights = fits.getdata(filename.replace(".fz", "_weight.fz"), 1)
            #dataerr = np.sqrt(1. / weights + data / gain)
            dataerr = 1. / weights + np.clip(data, 0, np.infty) / gain

            # Calculating flux density
            fnu = data * f0 * self.fnu_unit
            flam = fnu * const.c / wave**2
            flam = flam.to(self.flam_unit).value
            # Uncertaintis in flux density
            fnuerr = dataerr * f0 * self.fnu_unit
            flamerr = fnuerr * const.c / wave**2
            flamerr = flamerr.to(self.flam_unit).value
            flams.append(flam)
            flamerrs.append(flamerr)
            headers.append(h)
        flam = np.array(flams)
        flamerr = np.array(flamerrs)

        # Making new header with WCS
        wcs = WCS(h)
        add_before_ind = 2
        inds = [i + 1 for i in range(wcs.wcs.naxis)]
        inds.insert(add_before_ind, 0)
        newwcs = wcs.sub(inds)
        newwcs.wcs.ctype[add_before_ind] = ''
        newwcs.wcs.cname[add_before_ind] = ''
        # Making new header template
        newheader = headers[0].copy()
        newheader.update(newwcs.to_header())
        # Making table with metadata
        tab = []
        tab.append(self.bands)
        tab.append([self.wave_eff[band] for band in self.bands])
        names = ["FILTER", "WAVE_EFF"]
        hfields = ["OBJECT", "GAIN", "PSFFWHM", "DATE-OBS", "EXPTIME",
                   "EFECTIME", "NCOMBINE", "HIERARCH OAJ PRO FWHMMEAN"]
        for f in hfields:
            if not all([f in h for h in headers]):
                continue
            tab.append([h[f] for h in headers])
            names.append(f)
            if f in newheader:
                del newheader[f]
        tab = Table(tab, names=names)
        tab.rename_column("HIERARCH OAJ PRO FWHMMEAN", "PSFFWHM")
        tab.rename_column("OBJECT", "TILE")
        # Producing data cubes HDUs.
        hdus = [fits.PrimaryHDU()]
        hdu1 = fits.ImageHDU(flam, newheader)
        hdu1.header["EXTNAME"] = ("DATA", "Name of the extension")
        hdus.append(hdu1)
        hdu2 = fits.ImageHDU(flamerr, newheader)
        hdu2.header["EXTNAME"] = ("ERRORS", "Name of the extension")
        hdus.append(hdu2)
        for hdu in hdus:
            hdu.header["BSCALE"] = (1, "Linear factor in scaling equation")
            hdu.header["BZERO"] = (0, "Zero point in scaling equation")
            hdu.header["BUNIT"] = ("{}".format(self.flam_unit),
                                   "Physical units of the array values")
        thdu = fits.BinTableHDU(tab)
        hdus.append(thdu)
        thdu.header["EXTNAME"] = "METADATA"
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(self.cubename, overwrite=True)

    def make_cube(self, el=None, redo=False):
        if os.path.exists(self.cubename) and not redo:
            return
        cuts_exist = [os.path.exists(os.path.join(self.cutouts_dir, img)) for
                      img in self.cutnames]
        if not all(cuts_exist):
            print("Not all stamps are available locally, cube not produced.")
            return
        fields = [_.replace("_", "-") for _ in self.zps["FIELD"]]
        flams, flamerrs = [], []
        mask_star = []
        headers = []
        for band, img in zip(self.bands, self.cutnames):
            wave = self.wave_eff[band] * u.Angstrom
            filename = os.path.join(self.cutouts_dir, img)
            h = fits.getheader(filename, ext=1)
            tile = h["OBJECT"].replace("_", "-")
            # Getting zero-point calibration and location correction
            zp0 = self.zps[fields.index(tile)][band]
            x0 = h["X0TILE"]
            y0 = h["Y0TILE"]
            zpcorr = float(self.zpcorr[band](x0, y0))
            zp = zp0 + zpcorr
            gain = h["GAIN"]
            f0 = np.power(10, -0.4 * (48.6 + zp))
            data = fits.getdata(filename, 1)
            weights = fits.getdata(filename.replace(".fz", "_weight.fz"), 1)
            dataerr = np.sqrt(1 / weights + data / gain)

            # Calculating flux density
            fnu = data * f0 * self.fnu_unit
            flam = fnu * const.c / wave**2
            flam = flam.to(self.flam_unit).value
            # Uncertaintis in flux density
            fnuerr = dataerr * f0 * self.fnu_unit
            flamerr = fnuerr * const.c / wave**2
            flamerr = flamerr.to(self.flam_unit).value
            flams.append(flam)
            flamerrs.append(flamerr)
            headers.append(h)
        flam = np.array(flams)
        flamerr = np.array(flamerrs)

        # Making new header with WCS
        wcs = WCS(h)
        add_before_ind = 2
        inds = [i + 1 for i in range(wcs.wcs.naxis)]
        inds.insert(add_before_ind, 0)
        newwcs = wcs.sub(inds)
        newwcs.wcs.ctype[add_before_ind] = ''
        newwcs.wcs.cname[add_before_ind] = ''
        # Making new header template
        newheader = headers[0].copy()
        newheader.update(newwcs.to_header())
        # Making table with metadata
        tab = []
        tab.append(self.bands)
        tab.append([self.wave_eff[band] for band in self.bands])
        names = ["FILTER", "WAVE_EFF"]
        hfields = ["OBJECT", "GAIN", "PSFFWHM", "DATE-OBS", "EXPTIME",
                   "EFECTIME", "NCOMBINE", "HIERARCH OAJ PRO FWHMMEAN"]
        for f in hfields:
            if not all([f in h for h in headers]):
                continue
            tab.append([h[f] for h in headers])
            names.append(f)
            if f in newheader:
                del newheader[f]
        tab = Table(tab, names=names)
        tab.rename_column("HIERARCH OAJ PRO FWHMMEAN", "PSFFWHM")
        tab.rename_column("OBJECT", "TILE")
        # Producing data cubes HDUs.
        hdus = [fits.PrimaryHDU()]
        hdu1 = fits.ImageHDU(flam, newheader)
        hdu1.header["EXTNAME"] = ("DATA", "Name of the extension")
        hdus.append(hdu1)
        hdu2 = fits.ImageHDU(flamerr, newheader)
        hdu2.header["EXTNAME"] = ("ERRORS", "Name of the extension")
        hdus.append(hdu2)
        for hdu in hdus:
            hdu.header["BSCALE"] = (1, "Linear factor in scaling equation")
            hdu.header["BZERO"] = (0, "Zero point in scaling equation")
            hdu.header["BUNIT"] = ("{}".format(self.flam_unit),
                                   "Physical units of the array values")
        thdu = fits.BinTableHDU(tab)
        hdus.append(thdu)
        thdu.header["EXTNAME"] = "METADATA"
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(self.cubename, overwrite=True)

    def get_flam(self):
        if not os.path.exists(self.cubename):
            return None
        return fits.getdata(self.cubename, ext=1) * self.flam_unit

    def get_flamerr(self):
        if not os.path.exists(self.cubename):
            return None
        return fits.getdata(self.cubename, ext=2) * self.flam_unit

    def get_fnu(self):
        if not os.path.exists(self.cubename):
            return None
        flam = self.get_flam()
        fnu = (flam / const.c * self.wave[:, None, None]**2).to(
               self.fnu_unit)
        return fnu

    def get_fnuerr(self):
        if not os.path.exists(self.cubename):
            return None
        flamerr = self.get_flamerr()
        fnuerr = (flamerr / const.c * self.wave[:, None, None]**2).to(
                  self.fnu_unit)
        return fnuerr

    def get_mag(self):
        if not os.path.exists(self.cubename):
            return None
        fnu = self.get_fnu().value
        mab = -2.5 * np.log10(fnu) - 48.6
        return mab

    def get_magerr(self):
        if not os.path.exists(self.cubename):
            return None
        fnu = self.get_fnu().value
        fnuerr = self.get_fnuerr().value
        maberr = np.abs(2.5 / np.log(10) * fnuerr / fnu)
        return maberr

    def calc_el(self, el, output=None, store=True):
        """ Returns the EL fluxes """
        if not os.path.exists(self.cubename):
            return None
        flam = self.get_flam()
        flamerr = self.get_flamerr()

        if el=='halpha':
        # Calculating H-alpha
            wline = 6562.8 * u.AA
            el_bands = ["R", "F660", "I"]
            condition = [b in self.bands for b in el_bands]
            if not all(condition):
                print("Halpha estimation requires filters R, F660, I!")
                return
        if el=='oii':
        # Calculating [OII]
            wline = 3727 * u.AA
            el_bands = ["U", "F378", "G"]
            condition = [b in self.bands for b in el_bands]
            if not all(condition):
                print("[OII] estimation requires filters U, F378 and G!")
                return
        if self.survey=='S-PLUS':
            dir_filters = os.path.join(self.wdir,'S-PLUS_filters/filter_curves-master/')
            el_estimator = SPLUSEmLine(wline, el_bands, dir_filters)
        idx = np.array([self.bands.index(band) for band in el_bands])
        el_flux, el_flux_err = el_estimator(flam[idx], flamerr[idx])
        default_out = self.cubename.replace(".fits", f"_{el}.fits")
        output = default_out if output is None else output
        if store:
            # Saving fits
            h = fits.getheader(os.path.join(self.cutouts_dir, self.cutnames[0]),
                                            ext=1)
            h["EXTNAME"] = "DATA"
            hdu1 = fits.ImageHDU(el_flux.value, h)
            h["EXTNAME"] = "ERROR"
            hdu2 = fits.ImageHDU(el_flux_err.value, h)
            hdulist = fits.HDUList([fits.PrimaryHDU(), hdu1, hdu2])
            hdulist.writeto(output, overwrite=True)
        return el_flux, el_flux_err

    def dust_correction(self,halpha,el):
        wave_f660 = 6614.0 * u.Angstrom
        if not os.path.exists(self.cubename):
            return None
        fnu = self.get_fnu()
        mag_err = self.get_magerr()

        if el=='halpha':
            corr_bands = ["R", "F660", "I" ,"G"]
            condition = [b in self.bands for b in corr_bands]
            if not all(condition):
                print("Dust correction requires filters R, F660, I and G!")
        idx = np.array([self.bands.index(band) for band in corr_bands])
        fnus = fnu[idx].value
        magAB = -2.5 * np.log10(fnus) - 48.6
        g_i = magAB[3]-magAB[2]
        """" C = E(B-V) extinction law  with  eq. 20 or Vilella-Rojo+ (2015)"""
        ebv = np.array(0.206 * np.power(g_i, 1.68) - 0.0457)
        # ebv = excesso de cor
        ebv[np.isnan(ebv)] = 0
        x = 1 / wave_f660.to(u.micrometer).value
        wtran = (0.63 * u.micrometer)
        kappa = np.where(wave_f660 > wtran, 2.659 * (-1.857 + 1.040 * x),
                         2.659 * (-2.156 + 1.509 * x - 0.198 * x * x
                                  + 0.011 * (x * x * x)))
        return halpha * np.power(10, -0.4 * ebv * kappa)

    def nii_correction(self,halpha_nii, el):
        if not os.path.exists(self.cubename):
            return None
        mag = self.get_mag()
        if el=='halpha':
            corr_bands = ["I","G"]
            condition = [b in self.bands for b in corr_bands]
            if not all(condition):
                print("Dust correction requires filters G and I!")
        idx = np.array([self.bands.index(band) for band in corr_bands])
        magAB = mag[idx]
        #magAB_err = mag_err[idx]
        g_i = magAB[1]-magAB[0]

        """ Calculated NII corrected emission with eq. 21 or Vilella-Rojo+ (2015)"""
        idx = np.where(halpha_nii < 0)
        halpha = np.power(10, np.where((g_i <= 0.5),
                                       0.989 * np.log10(halpha_nii) - 0.193,
                                       0.954 * np.log10(halpha_nii) - 0.753))
        halpha[idx] = halpha_nii[idx]
        return halpha

    def color_diff(self, band1, band2):
        """ Returns color difference image """
        if not os.path.exists(self.cubename):
            return None
        flam = self.get_flam()
        flamerr = self.get_flamerr()
        condition = [b in self.bands for b in [band1, band2, "I"]]
        if not all(condition):
            print(f"Filter {band1} and {band2} required!")
            return
        idx = np.array([self.bands.index(band) for band in [band1, band2, "I"]])
        diff = (flam[idx[0]].value / flam[idx[1]].value)-1.
        diff2 = (flam[idx[0]].value / flam[idx[2]].value)-1.
        sel = (diff2<0.) #& (diff2<0.)
        error_band = flamerr[idx[0]].value
        diff[sel]=None
        plt.figure()
        cmap = plt.cm.get_cmap("gray")
        im = plt.imshow(diff, cmap=cmap, norm=LogNorm(), origin='lower')
        ylims = plt.gca().get_ylim()
        xlims = plt.gca().get_xlim()
        plt.ylim(ylims[0], ylims[1])
        plt.xlim(xlims[0], xlims[1])
        plt.axis('off')
        plt.tight_layout()
        plt.show()
        return



