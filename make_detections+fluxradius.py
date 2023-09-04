import sys
from splus_datacubes import SCube
import splusdata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sep
from astropy.io import fits
from matplotlib.patches import Ellipse, Circle
import matplotlib.cm as cm

###########  Input information:   ###############
# Directory where your image cutouts are going to be restore
dir_cutouts = '/your/cutout/directory/'
# Working directory, where filters and zeropoints are stored and results will be saved
dir_work = '/your/work/directory/'
# Survey of work (only S-PLUS for now)
survey = 'S-PLUS'
# If you are working with splusdata, add your credentials:
login_name =
password =
conn = splusdata.connect(login_name, password)

file_input = pd.read_csv('your-file-with-id-ra-dec-size.csv')
ids = file_input['Name']
ra = file_input['RA']
dec = file_input['DEC']
radius_img = file_input['radius']
fields = file_input['Field']

#Change to True if you wish to look at the R-image with the detected sources
show_extract_source= False

indexes = file_input.index.values.tolist()

deblend_var = np.array([0.05]*len(file_input))
deblend_var[0] = 0.005 ; deblend_var[8] = 0.005 ; deblend_var[9] = 0.005 ; deblend_var[10] = 0.005 ; deblend_var[19] = 0.005 ;
deblend_var[26] = 0.005 ; deblend_var[36] = 0.005 ; deblend_var[37] = 0.005 ; deblend_var[44] = 0.005 ; deblend_var[47] = 0.005;
deblend_var[48] = 0.005 ; deblend_var[51] = 0.005 ; deblend_var[57] = 0.005 ; deblend_var[59] = 0.005 ; deblend_var[60] = 0.005 ;
deblend_var[121] = 0.005 ; deblend_var[123] = 0.005 ; deblend_var[124] = 0.005 ; deblend_var[151] = 0.005 ; deblend_var[154] = 0.005
deblend_var[200] = 0.005 ; deblend_var[209] = 0.005 ; deblend_var[221] = 0.005; deblend_var[222] = 0.005; deblend_var[226] = 0.005

for i in indexes:
    # Galaxy information:
    galaxy = str(ids[i])
    print(i)
    print(galaxy)
    coords = [str(ra[i]), str(dec[i])]
    size = radius_img[i].item()
    field = fields[i]
    # Download stamps:
    scube = SCube(galaxy, coords, size, survey, conn=conn, field=field, wdir=dir_work,bands=['R','F660','I'],cut_dir=dir_cutouts)
    scube.download_stamps()

    data_r = fits.open(dir_cutouts + galaxy + '_R_' + str(size) + 'x' + str(size) + 'pix.fz')[1].data
    head_r = fits.open(dir_cutouts + galaxy + '_R_' + str(size) + 'x' + str(size) + 'pix.fz')[1].header

    bkg = sep.Background(data_r)
    bkg_image = bkg.back()

    data_sub = data_r - bkg
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms, deblend_cont=deblend_var[i])
    # To choose object with bigger a and closest to the center:
    a_array = objects['a']
    x_pos = objects['x']
    y_pos = objects['y']
    difference_array_x = np.absolute(x_pos - (size/2.))
    difference_array_y = np.absolute(y_pos - (size / 2.))

    #Try maximum a and check for its location:
    max_a_idx = a_array.argmax()
    idx_order = np.argsort(a_array)
    if (difference_array_x[max_a_idx]<25.) and (difference_array_y[max_a_idx]<25.):
        idx_gal = max_a_idx
    else:
        max_a_idx = idx_order[-2]
        if (difference_array_x[max_a_idx] < 25.) and (difference_array_y[max_a_idx] < 25.):
            idx_gal = np.argsort(a_array)[-2]
        else:
            idx_gal = np.argsort(a_array)[-3]

    idx_stars = np.array([True] * len(objects))
    idx_stars[idx_gal] = False

    stars_x = np.array(objects['x'][idx_stars])
    stars_y = np.array(objects['y'][idx_stars])
    stars_a = np.array(objects['a'][idx_stars])
    stars_b = np.array(objects['b'][idx_stars])
    stars_theta = np.array(objects['theta'][idx_stars])

    gal_a = np.array(objects['a'][~idx_stars])
    gal_x = np.array(objects['x'][~idx_stars])
    gal_y = np.array(objects['y'][~idx_stars])
    gal_b = np.array(objects['b'][~idx_stars])
    gal_theta = np.array(objects['theta'][~idx_stars])

    # Test object detection by plotting background-subtracted image with ellipses showing the detected objects
    # In blue, the galaxy of choice
    if show_extract_source == True:
        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, cmap='Greys_r', origin='lower', vmin=np.percentile(data_sub,0.1), vmax=np.percentile(data_sub,99.9))#, vmin=-0.1, vmax=3.5)

        # plot an ellipse for each object
        for i in range(len(stars_x)):
            e = Ellipse(xy=(stars_x[i], stars_y[i]),
                        width=6 * stars_a[i],
                        height=6 * stars_b[i],
                        angle=stars_theta[i] * 180. / np.pi)
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        cir = Circle(xy=(gal_x, gal_y), radius=4 * gal_a)
        cir.set_facecolor('none')
        cir.set_edgecolor('blue')
        ax.add_artist(cir)
        plt.show()
        pass

    else:
        star_mask = np.zeros(data_sub.shape, dtype=bool)
        sep.mask_ellipse(star_mask, stars_x, stars_y, 3 * stars_a, 3 * stars_b, stars_theta)
        print('Star mask done!')

        # Creating datacubes:
        scube.make_cube()
        print('Datacube done!')

        # Reading datacube:
        filename = dir_work + galaxy + '_' + str(size) + 'x' + str(
            size) + 'pix.fits'
        cube_data = fits.open(filename)
        flux_data = cube_data[1].data
        err_flux_data = cube_data[2].data

        flux_r = flux_data[0]
        flux_f660 = flux_data[1]
        flux_i = flux_data[2]
        fluxerr_r = err_flux_data[0]

        # Setting the circular radius that should contain all flux of the galaxy
        gal_radius_tot = 4 * gal_a[0]
        # Outer radius of the circular annulus to obtain the background estimation
        outer_anullus = gal_radius_tot + 3.

        data = flux_r.copy()
        data = data.byteswap(inplace=True).newbyteorder()
        dataerr = fluxerr_r.byteswap(inplace=True).newbyteorder()
        frac = np.array([0.2, 0.5, 0.7, 0.9])

        flux_tot, fluxerr_tot, flag = sep.sum_circle(data, gal_x, gal_y, gal_radius_tot, err=dataerr, mask=star_mask,
                                                     bkgann=(gal_radius_tot, outer_anullus))
        flux_radius, flags = sep.flux_radius(data, x=gal_x, y=gal_y, rmax=[gal_radius_tot], frac=frac, normflux=None,
                                             mask=star_mask)
        flux_radius = flux_radius[0]

        flux_r[star_mask] = None
        fig, ax = plt.subplots()
        cmap = cm.get_cmap('gray')
        im = ax.imshow(flux_r, interpolation='nearest', cmap=cmap, origin='lower', vmax=np.nanpercentile(flux_r, 99.1),
                       vmin=np.nanpercentile(flux_r, 0.1))
        ax.set_title(galaxy)
        ax.set_ylabel('DEC [pixel]')
        ax.set_xlabel('RA [pixel]')
        cmap.set_bad(color='black')
        cir = Circle(xy=(gal_x, gal_y), radius=gal_radius_tot)
        cir.set_facecolor('none')
        cir.set_edgecolor('blue')
        cir.set_label('Total Flux')
        ax.add_artist(cir)
        cir0 = Circle(xy=(gal_x, gal_y), radius=outer_anullus)
        cir0.set_facecolor('none')
        cir0.set_edgecolor('blue')
        cir0.set_linestyle('--')
        ax.add_artist(cir0)
        colors = ['green', 'magenta', 'cyan', 'purple']
        labels = ['20% Flux', '50% Flux', '70% Flux', '90% Flux']
        for f in range(0, len(flux_radius)):
            cir1 = Circle(xy=(gal_x, gal_y), radius=flux_radius[f])
            cir1.set_facecolor('none')
            cir1.set_edgecolor(colors[f])
            cir1.set_label(labels[f])
            ax.add_artist(cir1)
        plt.legend()
        plt.savefig(dir_work +'/jelly_radius_figs/{}_fluxradius.png'.format(galaxy))

        file_input['radius_20'][i] = flux_radius[0]
        file_input['radius_50'][i] = flux_radius[1]
        file_input['radius_70'][i] = flux_radius[2]
        file_input['radius_90'][i] = flux_radius[3]

        file_input.to_csv('/your-file-with-id-ra-dec-size.csv',index=False)
        print('File updated with flux radius!')

