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
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)

    # To choose object with bigger a and closest to the center:
    a_array = objects['a']
    x_pos = objects['x']
    y_pos = objects['y']
    difference_array_x = np.absolute(x_pos - (size/2.))
    difference_array_y = np.absolute(y_pos - (size / 2.))

    #Try maximum a and check for its location:
    max_a_idx = a_array.argmax()
    if (difference_array_x[max_a_idx]<12.) and (difference_array_y[max_a_idx]<12.):
        idx_gal = max_a_idx
    else:
        idx_gal = np.argsort(a_array)[-2]

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
        im = ax.imshow(data_sub, cmap='Greys_r', origin='lower', vmin=-0.1, vmax=3.5)

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