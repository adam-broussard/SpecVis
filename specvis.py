import numpy as np
from pandas import read_csv
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Cursor
import matplotlib.patheffects as PathEffects
import fsps
from glob import glob
import photosim

phot = photosim.photosim()


def plot_spectra(spread = 'dust', tau = 0.1, Av = 0.0, plotfilters = [72, 73, 74, 75, 76, 160, 161, 162], filtnames = ['u', 'g', 'r', 'i', 'z', 'J', 'H', 'K']):

    fig, sp = plt.subplots(figsize = (8,8))
    plt.subplots_adjust(left=0.25, bottom=0.25)

    clone = sp.twinx()

    names, waves, throughputs = read_filters()

    names = [names[x] for x in plotfilters]
    waves = [waves[x] for x in plotfilters]
    throughputs = [throughputs[x] for x in plotfilters]

    for x, (thisname, filtwave, filttrans) in enumerate(zip(filtnames, waves, throughputs)):

        sp.plot(filtwave, filttrans, color = plt.cm.RdYlBu_r(int(float(x)*plt.cm.RdYlBu.N/len(waves))), label = thisname)
        txt = sp.text(np.average(filtwave), max(filttrans) + .1, thisname, color = plt.cm.RdYlBu_r(int(float(x)*plt.cm.RdYlBu.N/len(waves))), fontsize = 16, ha = 'center')
        txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])



    zax = plt.axes([0.25, 0.1, 0.65, 0.03])

    RYB = plt.get_cmap('RdYlBu')

    lines = []

    if spread == 'dust':

        dustvals = np.linspace(0, 5, 6)

        for x, thisdust in enumerate(dustvals):
            wave, spec = phot.find_spectrum(tage = 2., sfh_params = {'tau': tau}, sfh_type = 6, peraa = True, Av = thisdust)
            lines.append(clone.plot(wave, np.log10(spec), color = plt.cm.RdYlBu_r(int(float(x)*plt.cm.RdYlBu_r.N/len(dustvals))))[0])

        clone.set_ylim(40, 50)

    elif spread == 'sfh' or spread == 'sfr':

        tauvals = np.append(-np.logspace(-1,0,3), np.logspace(-1,0,3)[::-1])

        for x, thistau in enumerate(tauvals):
            wave, spec = phot.find_spectrum(tage = 2., sfh_params = {'tau': thistau}, sfh_type = 6, peraa = True, Av = Av)
            lines.append(clone.plot(wave, np.log10(spec), color = plt.cm.RdYlBu(int(float(x)*plt.cm.RdYlBu.N/len(tauvals))))[0])

        clone.set_ylim(30, 50)

    zslider = Slider(zax, 'Redshift', 0.0, 3.0, valinit=0.0)

    origx = [thisline.get_xdata() for thisline in lines]

    def update(val):
        redshift = zslider.val
        for thisline, thisoriginal_x in zip(lines, origx):
            thisline.set_xdata(thisoriginal_x * (1.+redshift))
        fig.canvas.draw_idle()
    zslider.on_changed(update)

    sp.set_xlabel(r'Wavelength ($\AA$)')
    sp.set_ylabel('Filter Throughput')
    sp.set_xscale('log')
    sp.set_xlim(10**3, 10**4.5)
    sp.set_ylim(0, 5)
    clone.set_ylabel('log(Flux)')

    return zslider



def plot_spectra_interactive(plotfilters = [72, 73, 74, 75, 76, 160, 161, 162], filtnames = ['u', 'g', 'r', 'i', 'z', 'J', 'H', 'K']):

    fig, sp = plt.subplots(figsize = (10,10))
    plt.subplots_adjust(left=0.15, bottom=0.3, right = 0.85)

    clone = sp.twinx()

    names, waves, throughputs = read_filters()

    names = [names[x] for x in plotfilters]
    waves = [waves[x] for x in plotfilters]
    throughputs = [throughputs[x] for x in plotfilters]

    for x, (thisname, filtwave, filttrans) in enumerate(zip(filtnames, waves, throughputs)):

        sp.plot(filtwave, filttrans, color = 'k', linewidth = 3.5)
        sp.plot(filtwave, filttrans, color = plt.cm.RdYlBu_r(int(float(x)*plt.cm.RdYlBu.N/len(waves))), label = thisname)
        txt = sp.text(np.average(filtwave), max(filttrans) + .1, thisname, color = plt.cm.RdYlBu_r(int(float(x)*plt.cm.RdYlBu.N/len(waves))), fontsize = 16, ha = 'center')
        txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])



    zax = plt.axes([0.1, 0.15, 0.8, 0.03])
    Avax = plt.axes([0.1, 0.1, 0.8, 0.03])
    tauax = plt.axes([0.1, 0.05, 0.8, 0.03])

    zslider = Slider(zax, 'Redshift', 0.0, 3.0, valinit=0.0)
    Avslider = Slider(Avax, '$A_v$', 0.0, 5., valinit=0.0)
    tauslider = Slider(tauax, r'$\tau$', -1, 1, valinit=-1, valstep = 0.02)

    wave, spec = phot.find_spectrum(tage = 2., sfh_params = {'tau': -1}, sfh_type = 6, peraa = True, Av = 0.0)
    line = clone.plot(wave, np.log10(spec), color = 'k', zorder = 1)[0]

    photo_x, photo_y = get_photometry(wave, spec, waves, throughputs)
    photometry = clone.scatter(photo_x, np.log10(photo_y), s = 75, edgecolor = 'k', linewidth = 0.5, facecolor = plt.cm.RdYlBu_r(np.linspace(0, plt.cm.RdYlBu_r.N, len(photo_x), dtype = int)), zorder = 2)

    def update(val):
        redshift = zslider.val
        tau = tauslider.val
        Av = Avslider.val
    
        newx, newy = phot.find_spectrum(tage = 2., sfh_params = {'tau': tau}, sfh_type = 6, peraa = True, Av = Av)
        newx = newx * (1.+redshift)
        line.set_xdata(newx)
        line.set_ydata(np.log10(newy))

        newphoto_x, newphoto_y = get_photometry(newx, newy, waves, throughputs)
        photometry.set_offsets(list(zip(newphoto_x, np.log10(newphoto_y))))

        avgy = np.average(np.log10(newphoto_y))
        clone.set_ylim(avgy-10, avgy+10)

        fig.canvas.draw_idle()

    zslider.on_changed(update)
    Avslider.on_changed(update)
    tauslider.on_changed(update)

    sp.set_xlabel(r'Wavelength ($\AA$)')
    sp.set_ylabel('Filter Throughput')
    sp.set_xscale('log')
    sp.set_xlim(10**3, 10**5)
    sp.set_ylim(0, 5)
    clone.set_ylabel('log(Flux)')
    clone.set_ylim(30, 50)

    return zslider, Avslider, tauslider




def read_filters(fpath = './cosmos_example.FILTER.RES.DR34.latest'):

    names = []
    waves = []
    throughputs = []

    with open(fpath, 'r') as readfile:
        all_lines = readfile.readlines()

    line_num = 0

    while line_num < len(all_lines):
        linesplit = all_lines[line_num].split()
        filt_len = int(linesplit[0])
        names.append(linesplit[1])
        tempwave, tempthrough = np.array(list(zip(*[thisline[:-1].split(' ')[-2:] for thisline in all_lines[line_num+1:line_num+1+filt_len]])), dtype = float)
        throughputs.append(tempthrough)
        waves.append(tempwave)
        line_num += filt_len + 1

    # for x in range(len(waves)):

    #     throughputs[x] = throughputs[x]/np.trapz(throughputs[x], x = waves[x])

    return names, waves, throughputs



def get_photometry(wavelength, f_lambda, filterwaves, filterthroughs, wavestep = 0.5, f_nu = False):

    c = 3e10 # Speed of Light in cgs

    interp_x = np.arange(3000, 24000, wavestep) # The x-axis of the spectrum

    phot = []
    phot_wave = []

    for x in range(len(filterwaves)):
        filterthroughs[x] = filterthroughs[x]/np.trapz(filterthroughs[x], x = filterwaves[x])

    # Loop through each of the LSST filters
    for index, (thiswave, thisthrough) in enumerate(zip(filterwaves, filterthroughs)):

        interp_y = np.interp(interp_x, wavelength, f_lambda) # The interpolated y-axis of the spectrum (in f_lambda)

        # Find the filter curve values at the same x values as interp_x
        filter_interp_y = np.interp(interp_x, thiswave, thisthrough) 

        phot.append(10**-8 * np.trapz(filter_interp_y * interp_y * interp_x, x = interp_x)/np.trapz(filter_interp_y * c / interp_x, x = interp_x))
        phot_wave.append(np.trapz(thisthrough * thiswave, x = thiswave))

    phot_wave = np.array(phot_wave)
    phot = np.array(phot)

    if f_nu:
        return phot_wave, phot
    else:
        return phot_wave, phot * c/phot_wave**2 * 10**8