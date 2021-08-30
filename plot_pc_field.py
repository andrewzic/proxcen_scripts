import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation
from astropy.wcs import WCS
from astropy.time import Time
from astropy.nddata.utils import Cutout2D

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle

import numpy as np


p = SkyCoord(217.38469608*u.deg, -62.67536234*u.deg) #pm-corrected coords
#217.39346574260*u.deg, -62.67618210292*u.deg) #old coords

# t = Time('2019-05-02 08:00:00', format = 'iso')

# p_new,res = astrometry(p, t)
# print(p_new.to_string('hmsdms', sep = ':'))
# print(p_new)

# proxima_coord = p_new


field_imagef = 'proxima_sb8604_field_I.image.tt0.fits'
before_imagef = 't_0007_I.image.fits'
flare_imagef = 't_0009_I.image.fits'

hdu_field = fits.open(field_imagef)[0]
wcs_field = WCS(hdu_field.header, naxis = 2)

hdu_before = fits.open(before_imagef)[0]
hdu_flare = fits.open(flare_imagef)[0]
wcs_flare = WCS(hdu_flare.header, naxis = 2)
wcs_before = WCS(hdu_before.header, naxis = 2)

fig = plt.figure(figsize = (9/1.6, 12/1.6))#, wspace = 0)#, hspace = 0)

ax_field = plt.subplot2grid((3, 2), (0, 0), rowspan = 2, colspan = 2, fig = fig, projection = wcs_field)
ax_before = plt.subplot2grid((3,2), (2, 0), fig = fig, projection = wcs_before)
ax_flare = plt.subplot2grid((3,2), (2, 1), fig = fig, projection = wcs_flare)#, sharex = ax_before)#

ax_before.yaxis.set_label_position('right')
ax_flare.yaxis.set_label_position('right')

ax_before.yaxis.tick_right()
ax_flare.yaxis.tick_right()



plt.sca(ax_field)
print(hdu_field.data.shape)
ax_field_im = ax_field.imshow(1E6*hdu_field.data[0][0], vmin = -0.00015*1E6, vmax = 0.0004*1E6, cmap = 'gray_r', origin = 'lower')
ax_field.set_xlabel('Right Ascension (J2000)', fontsize = 14)
ax_field.set_ylabel('Declination (J2000)', fontsize = 14)


boxwidth = 20.0 #width of zoom box in arcminutes

data_before = Cutout2D(hdu_before.data[0][0], p, (35.0*u.arcminute, 35.0*u.arcminute), wcs = wcs_before).data
data_flare = Cutout2D(hdu_flare.data[0][0], p, (35.0*u.arcminute, 35.0*u.arcminute), wcs = wcs_flare).data

middle_ind = data_before.shape[0]//2

ax_before_im = ax_before.imshow(1E3*hdu_before.data[0][0], vmin = -0.02*1E3, vmax = 0.11*1E3, cmap = 'gray_r', aspect = 'auto', origin = 'lower')
print('test')

ax_flare_im = ax_flare.imshow(1E3*hdu_flare.data[0][0], vmin = -0.02*1E3, vmax = 0.11*1E3, cmap = 'gray_r', aspect = 'auto', origin = 'lower')

top_left = np.array([[(p.ra + boxwidth/np.cos(p.dec)/2.0*u.arcminute).value, (p.dec - 1.0*boxwidth/2.0*u.arcminute).value]])
top_left_pix = wcs_before.wcs_world2pix(top_left, 0)
bottom_right = np.array([[(p.ra - boxwidth/np.cos(p.dec)/2.0*u.arcminute).value, (p.dec + 1.0*boxwidth/2.0*u.arcminute).value]])
bottom_right_pix = wcs_before.wcs_world2pix(bottom_right, 0)


ax_before.set_xlim(top_left_pix[0,0], bottom_right_pix[0,0])#, transform = ax_field.get_transform('world'))
ax_before.set_ylim( top_left_pix[0,1], bottom_right_pix[0,1]+1.0)


top_left = np.array([[(p.ra + boxwidth/np.cos(p.dec)/2.0*u.arcminute).value, (p.dec - 1.0*boxwidth/2.0*u.arcminute).value]])
top_left_pix = wcs_flare.wcs_world2pix(top_left, 0)
bottom_right = np.array([[(p.ra - boxwidth/np.cos(p.dec)/2.0*u.arcminute).value, (p.dec + 1.0*boxwidth/2.0*u.arcminute).value]])
bottom_right_pix = wcs_flare.wcs_world2pix(bottom_right, 0)


ax_flare.set_xlim(top_left_pix[0,0], bottom_right_pix[0,0])#, transform = ax_field.get_transform('world'))
ax_flare.set_ylim( top_left_pix[0,1], bottom_right_pix[0,1]+1.0)

#ax_before.set(xlim = ((p.ra - 0.5*u.arcminute).value, (p.ra + 0.5*u.arcminute).value),
#              ylim = ((p.dec - 0.5*u.arcminute).value, (p.dec + 0.5*u.arcminute).value))
#ax_flare.set(xlim = ((p.ra - 0.5*u.arcminute).value, (p.ra + 0.5*u.arcminute).value),
#              ylim = ((p.dec - 0.5*u.arcminute).value, (p.dec + 0.5*u.arcminute).value))




for ax in fig.get_axes():
    ax.set_aspect('equal', 'box')

#plt.tight_layout()
ax_before.set_xticklabels([])

bottom_left = np.array([[(p.ra + boxwidth/np.cos(p.dec)/2.0*u.arcminute).value, (p.dec - boxwidth/2.0*u.arcminute).value]])
bottom_left_pix = wcs_field.wcs_world2pix(bottom_left, 0)[0]

top_right = np.array([[(p.ra - boxwidth/np.cos(p.dec)/2.0*u.arcminute).value, (p.dec + boxwidth/2.0*u.arcminute).value]])
top_right_pix = wcs_field.wcs_world2pix(top_right, 0)[0]

w,h = top_right_pix - bottom_left_pix
print(w,h)

r = Rectangle(bottom_left_pix, w, h, edgecolor = 'red', facecolor = 'None', linewidth = 1.5)
ax_field.add_patch(r)


ra, dec = ax_field.coords
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_ticks_position('tb')
ra.set_ticklabel_position('t')
ra.set_axislabel_position('t')

ax_field.tick_params('both', direction = 'in', color = 'black')
ax_field.tick_params('both', which = 'major', length = 6)
ax_field.tick_params('both', which = 'minor', length = 3)

mid_pix = np.array([[p.ra.value, p.dec.value]])
mid_pix = wcs_before.wcs_world2pix(mid_pix, 0)[0]
#ax_field.scatter([p.ra.value], [p.dec.value], marker = 'o', facecolor = 'white', s = 25, transform = ax_field.get_transform('world'))


wcss = [wcs_flare, wcs_before]
for i_, ax in enumerate([ax_flare, ax_before]):
    plt.sca(ax)
    ra_offs, dec_offs = ax.get_coords_overlay(p.skyoffset_frame())
    plt.minorticks_on()
    ra_offs.set_coord_type('longitude', 180)
    ra_offs.set_format_unit(u.arcminute, decimal = True)
    dec_offs.set_format_unit(u.arcminute, decimal = True)
    ra_offs.tick_params(direction = 'in', color = 'black')
    ra_offs.tick_params(which = 'major', length = 6)
    ra_offs.tick_params(which = 'minor', length = 3)

    dec_offs.tick_params(direction = 'in', color = 'black')
    dec_offs.tick_params(which = 'major', length = 6)
    dec_offs.tick_params(which = 'minor', length = 3)

    #ra_offs.set_ticklabel_position('b')


    ax.plot([mid_pix[0], mid_pix[0]], mid_pix[1]*np.array([1.02, 1.06]), c = 'red', linewidth = 2.0)#, transform = ax.get_transform('world'))
    ax.plot(mid_pix[0]*np.array([1.02, 1.06]), [mid_pix[1], mid_pix[1]], c = 'red', linewidth = 2.0)#, transform = ax.get_transform('world'))
    ra, dec = ax.coords

    ra.set_ticks_visible(False)
    ra.set_ticklabel_visible(False)
    dec.set_ticks_visible(False)
    dec.set_ticklabel_visible(False)
    ra_offs.display_minor_ticks(True)
    dec_offs.display_minor_ticks(True)
    if i_ == 0:
        #dec.set_ticks_visible(False)
        dec_offs.set_ticklabel_visible(False)
    else:
        dec_offs.set_ticks_position('lr')
        dec_offs.set_ticklabel_position('l')
        dec_offs.set_axislabel_position('l')
        dec_offs.set_axislabel('Dec. Offset', fontsize = 14)

#         #ra.set_ticks_visible(False)
#         ra.set_ticklabel_visible(False)
#         ra.set_axislabel('')
#     else:
#         ra.set_ticklabel_position('b')
#         ra.set_axislabel_position('b')

    ra_offs.set_ticks_position('tb')
    ra_offs.set_ticklabel_position('b')
    ra_offs.set_axislabel_position('b')
    ra_offs.set_axislabel('R.A. Offset', fontsize = 14)

#     ax.tick_params('both', direction = 'in', color = 'white')
#     ax.tick_params('both', which = 'major', length = 5)
#     ax.tick_params('both', which = 'minor', length = 3)

    boxwidth = 4.0
    wcs_ = wcss[i_]
    bottom_left = np.array([[(p.ra + boxwidth/2.0*u.arcminute).value, (p.dec - np.cos(p.dec)*boxwidth/2.0*u.arcminute).value]])
    bottom_left_pix = wcs_.wcs_world2pix(bottom_left, 0)[0]

    top_right = np.array([[(p.ra - boxwidth/2.0*u.arcminute).value, (p.dec + np.cos(p.dec)*boxwidth/2.0*u.arcminute).value]])
    top_right_pix = wcs_.wcs_world2pix(top_right, 0)[0]

    w,h = top_right_pix - bottom_left_pix
    print(w,h)

    r = Rectangle(bottom_left_pix, w, h, edgecolor = 'red', facecolor = 'None', linewidth = 0.75, linestyle = '--')
#     ax.add_patch(r)

# ax_before.set_ylabel('Declination (J2000)', fontsize = 14)
# ax_before.set_xlabel('Right Ascension (J2000)', fontsize = 14)
# ax_flare.set_xlabel('Right Ascension (J2000)', fontsize = 14)

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.0, hspace=0.0)

ax_field_pos = ax_field.get_position()

cb_ax = fig.add_axes([ax_field_pos.xmax + 0.00, ax_field_pos.ymin, 0.03, ax_field_pos.ymax - ax_field_pos.ymin])
cbar = fig.colorbar(ax_field_im, cax=cb_ax)
cbar.set_label(r'Flux Density ($\mu$Jy$~$beam$^{-1}$)', fontsize = 12)

ax_flare_pos = ax_flare.get_position()

cb_ax = fig.add_axes([ax_flare_pos.xmax + 0.00, ax_flare_pos.ymin, 0.03, ax_flare_pos.ymax - ax_flare_pos.ymin])
cbar = fig.colorbar(ax_flare_im, cax=cb_ax)
cbar.set_label(r'Flux Density (mJy$~$beam$^{-1}$)', fontsize = 12)

t_before = Time(hdu_before.header['DATE-OBS'], format = 'isot', location = EarthLocation.of_site('MRO'))
t_flare = Time(hdu_flare.header['DATE-OBS'], format = 'isot', location = EarthLocation.of_site('MRO'))
times = Time([t_before, t_flare], location = EarthLocation.of_site('MRO'))


ltt_before = t_before.light_travel_time(p)
ltt_bary = np.array([t.light_travel_time(p) for t in times])

times = times.tdb + ltt_bary

t_before = times[0]
t_flare = times[1]

ax_flare.text(0.08, 0.86, '{:.5f} MBJD'.format(t_flare.mjd), fontsize = 12, bbox=dict(facecolor='white', alpha=0.5), transform=ax_flare.transAxes)

ax_before.text(0.08, 0.86, '{:.5f} MBJD'.format(t_before.mjd), fontsize = 12, bbox=dict(facecolor='white', alpha=0.5), transform=ax_before.transAxes)

#plt.tight_layout()


#fig.show_markers(ra, dec, layer = 'markers', edgecolor = 'red', facecolor = 'none', s = 3000, marker = path, linewidth = 1.8)
plt.savefig('proxima_field_flare_map_gr.pdf', bbox_inches = 'tight')
plt.savefig('proxima_field_flare_map_gr.png', bbox_inches = 'tight', dpi = 300, facecolor = 'w')
plt.show()
