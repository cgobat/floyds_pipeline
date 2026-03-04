
import datetime
import glob
import os
import re
import string
import sys
import time

import numpy as np
from pyraf import iraf
from astropy.io import fits

import floyds

pyversion = sys.version_info[0]

def sortbyJD(lista):

    JDlist = []
    for img in lista:
        hdr = fits.open(img)[0].header
        if 'MJD' in hdr:
            JDlist.append(hdr['MJD'])
        else:
            JDlist.append(hdr['MJD-OBS'])
    lista = np.array(lista)
    JDlist = np.array(JDlist)
    inds = JDlist.argsort()
    sortedlista = lista[inds]
    return list(sortedlista)

def ask(question):
    if pyversion>=3:
        answ = input(question)
    else:
        answ = raw_input(question)
    return answ


# ###########################################################
def ReadAscii2(ascifile):

    f = open(ascifile, 'r')
    ss = f.readlines()
    f.close()
    vec1, vec2 = [], []
    for line in ss:
        if line[0] != '#':
            vec1.append(float(str.split(line)[0]))
            vec2.append(float(str.split(line)[1]))
    return vec1, vec2


#########################################################################
def readspectrum(img):

    fl = ''
    lam = ''
    graf = 1
    spec = fits.open(img)
    head = spec[0].header
    try:
        if spec[0].data.ndim == 1:
            fl = spec[0].data
        elif spec[0].data.ndim == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.ndim == 3:
            fl = spec[0].data[0, 0, :]
    except:
        if spec[0].data.rank == 1:
            fl = spec[0].data
        elif spec[0].data.rank == 2:
            fl = spec[0].data[:, 0]
        elif spec[0].data.rank == 3:
            fl = spec[0].data[0, 0, :]
    naxis1 = head['naxis1']
    try:
        crpix1 = head['crpix1']
        crval1 = head['crval1']
        try:
            cdelt1 = head['cdelt1']
        except:
            cdelt1 = head['cd1_1']
        pix = np.array(range(1, naxis1 + 1, 1))
        pix = np.array(range(1, len(fl) + 1, 1))
        lam = (pix - crpix1) * cdelt1 + crval1
    except:
        try:
            WAT = head['WAT2_001']
            pix = np.array(range(1, naxis1 + 1, 1))
            crpix1 = str.split(str.split(WAT, '"')[1])[0]
            crval1 = str.split(str.split(WAT, '"')[1])[3]
            cdelt1 = str.split(str.split(WAT, '"')[1])[4]
            lam = (pix - float(crpix1)) * float(cdelt1) + float(crval1)
        except:
            graf = 0
    return lam, fl


###########################################################################

def readlist(listfile):

    if '*' in listfile:
        imglist = glob.glob(listfile)
    elif ',' in listfile:
        imglist = str.split(listfile, ',')
    else:
        try:
            hdulist = fits.open(listfile)
        except:
            hdulist = []
        if hdulist:
            imglist = [listfile]
        else:
            try:
                ff = open(listfile, 'r')
                files = ff.readlines()
                ff.close()
                imglist = []
                for ff in files:
                    ff = re.sub(' ', '', ff)
                    if not ff == '\n' and ff[0] != '#':
                        ff = re.sub('\n', '', ff)
                        try:
                            hdulist = fits.open(ff)
                            imglist.append(ff)
                        except Exception as e:
                            print('problem reading header of', ff)
                            print(e)
            except:
                sys.exit('\n##### Error ###\n file ' + str(listfile) + ' do not  exist\n')
    if len(imglist) == 0:
        sys.exit('\n##### Error ###\nIf "' + str(listfile) \
                 + '" is an image, it is corrupted \n or is not a list of image\n')
    return imglist


##############################################################################
def delete(listfile):

    if listfile[0] == '@':
        ff = open(listfile[1:])
        files = ff.readlines()
        imglist = []
        for ff in files:
            ff = re.sub(' ', '', ff)
            if not ff == '\n' and ff[0] != '#':
                ff = re.sub('\n', '', ff)
                imglist.append(ff)
    elif ',' in listfile:
        imglist = str.split(listfile, ',')
    else:
        imglist = [listfile]
    lista = []
    for _file in imglist:   lista = lista + glob.glob(_file)
    if lista:
        for _file in lista:
            if os.path.isfile(_file):
                try:
                    os.system('rm ' + _file)
                except:
                    pass


###############################################################
def readhdr(img):
    try:
        hdr = fits.getheader(img)
    except Exception as e:
        print("Couldn't read header of {}. Try deleting it and starting over.".format(img))
        raise e
    return hdr

def readkey3(hdr, keyword):

    try:
        _instrume = hdr.get('INSTRUME').lower()
    except:
        _instrume = 'none'
    useful_keys = {'object': 'OBJECT',
                       'date-obs': 'DATE-OBS'}
    if _instrume in ['en05', 'en06', 'en12']:
        useful_keys.update(
            {'ut': 'UTSTART',
             'obstype': 'OBSTYPE',
             'RA': 'RA',
             'DEC': 'DEC',
             'CAT-RA': 'CAT-RA',
             'CAT-DEC': 'CAT-DEC',
             'datamin': -100,
             'datamax': 60000,
             'grpid': 'BLKUID',
             'exptime': 'EXPTIME'
             })
        if not hdr.get('HDRVER'):
            useful_keys.update(
                {'JD': 'MJD',
                 'MJD': 'MJD',
                 'lamp': 'LMP_ID',
                 'gain': 'GAIN',
                 'instrume': 'INSTRUME',
                 'grism': 'GRISM',
                 'ron': 'READNOIS',
                 'airmass': 'AIRMASS',
                 'slit': 'APERWID',
                 'telescop': 'TELESCOP'
                 })
        else:
            useful_keys.update(
                {'JD': 'MJD-OBS',
                 'MJD': 'MJD-OBS',
                 'lamp': 'LMP1ID',
                 'gain': 'GAIN',
                 'instrume': 'INSTRUME',
                 'grism': 'GRISM',
                 'ron': 'RDNOISE',
                 'airmass': 'AIRMASS',
                 'slit': 'APERWID',
                 'telescop': 'TELESCOP'
                 })

    if keyword in useful_keys:
        if type(useful_keys[keyword]) == float:
            value = useful_keys[keyword]
        else:
            value = hdr.get(useful_keys[keyword])
            if keyword == 'date-obs':

                try:
                    value = re.sub('-', '', str.split(value, 'T')[0])
                except:
                    pass
            elif keyword == 'ut':

                try:
                    value = str.split(value, 'T')[1]
                except:
                    pass
            elif keyword == 'JD':
                value = float(value) + 0.5
            elif keyword == 'instrume':
                value = value.lower()
            elif keyword == 'grism':
                if not value: value = 'full'
            elif keyword == 'grpid':
                value = '_'.join([hdr.get('BLKUID'), hdr.get('OBJECT')])
            elif keyword == 'RA' or keyword == 'CAT-RA':

                value0 = str.split(value, ':')
                try:
                    value = ((float(value0[0]) + ((float(value0[1]) + (float(value0[2]) / 60.)) / 60.)) * 15)
                except:
                    value = 0
            elif keyword == 'DEC' or keyword == 'CAT-DEC':

                value0 = str.split(value, ':')
                try:
                    if '-' in str(value0[0]):
                        value = ((-1) * (np.abs(float(value0[0])) + ((float(value0[1]) + (float(value0[2]) / 60.)) / 60.)))
                    else:
                        value = (float(value0[0]) + ((float(value0[1]) + (float(value0[2]) / 60.)) / 60.))
                except:
                    value = 0
            elif keyword == 'slit':
                value = re.sub('\"', '', re.sub('slit', '', str(value).lower()))
                value = re.sub('as', '', re.sub('_', '', str(value.lower())))
                if value == 'UNKNOWN': value = '1.6'
            elif keyword == 'object':
                value = re.sub('\}', '', value)
                #value = re.sub('\+', '', value)
                value = re.sub('\{', '', value)
                value = re.sub('\[', '', value)
                value = re.sub('\]', '', value)
                value = re.sub('\(', '', value)
                value = re.sub('\)', '', value)
                value = re.sub('\,', '', value)
                value = re.sub(' ', '', value)
    else:
        if keyword == 'date-night':

            _date = readkey3(hdr, 'DATE-OBS')
            a = (datetime.datetime.strptime(str.split(_date, '.')[0], "20%y-%m-%dT%H:%M:%S") - datetime.timedelta(
                .0)).isoformat()
            value = re.sub('-', '', str.split(a, 'T')[0])
        elif keyword == 'TELID':
            value = hdr.get(keyword)
            value = re.sub('-', '', value)
            if value not in ['fts', 'ftn']:
                value = hdr.get('SITEID')
                if value in ['ogg']: value = 'ftn'
                if value in ['coj']: value = 'fts'
        elif keyword == 'PROPID':
            value = hdr.get('PROPID')
            # Remove whitespace from proposal id
            value = re.sub(r"\s+", "", value)
        else:
            try:
                value = hdr.get(keyword)
            except:
                sys.exit('Warning: keyword not valid')
    if type(value) == str:    value = re.sub('\#', '', value)
    return value


#######################################################
def writeinthelog(text, logfile):
    f = open(logfile, 'a')
    f.write(text)
    f.close()

def updateheader(filename, dimension, headerdict):
    tupledict = {key: tuple(value) for key, value in headerdict.items()}
    try:
        hdulist = fits.open(filename, mode='update')
        header = hdulist[dimension].header
        header.update(tupledict)
        hdulist.close()
    except Exception as e:
        print('header of', image, 'not updated:')
        print(e)

#################################################################################################
def display_image(img, frame, _z1, _z2, scale, _xcen=0.5, _ycen=0.5, _xsize=1, _ysize=1, _erase='yes'):
    goon = 'True'

    ds9 = subprocess.Popen("ps -U" + str(os.getuid()) + "|grep -v grep | grep ds9", shell=True,
                           stdout=subprocess.PIPE).stdout.readlines()
    if len(ds9) == 0:
        subproc = subprocess.Popen('ds9', shell=True)
        time.sleep(3)

    if glob.glob(img):

        iraf.images(_doprint=0)
        iraf.tv(_doprint=0)

        if _z2:
            try:
                sss = iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,
                                   fill='yes', zscale='no', zrange='no', z1=_z1, z2=_z2, Stdout=1)
            except:
                print('')
                print('### ERROR: PROBLEM OPENING DS9')
                print('')
                goon = 'False'
        else:
            try:
                sss = iraf.display(img, frame, xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize, erase=_erase,
                                   fill='yes', Stdout=1)
            except:
                print('')
                print('### ERROR: PROBLEM OPENING DS9')
                print('')
                goon = False

        if scale and goon:
            answ0 = ask('>>> Cuts OK ? [y/n] ? [y] ')
            if not answ0:
                answ0 = 'y'
            elif answ0 == 'no' or answ0 == 'NO':
                answ0 = 'n'

            while answ0 == 'n':
                _z11 = float(str.split(str.split(sss[0])[0], '=')[1])
                _z22 = float(str.split(str.split(sss[0])[1], '=')[1])
                z11 = ask('>>> z1 = ? [' + str(_z11) + '] ? ')
                z22 = ask('>>> z2 = ? [' + str(_z22) + '] ? ')
                if not z11:
                    z11 = _z11
                else:
                    z11 = float(z11)
                if not z22:
                    z22 = _z22
                else:
                    z22 = float(z22)
                print(z11, z22)
                sss = iraf.display(img, frame, fill='yes', xcen=_xcen, ycen=_ycen, xsize=_xsize, ysize=_ysize,
                                   erase=_erase,
                                   zrange='no', zscale='no', z1=z11, z2=z22, Stdout=1)
                answ0 = ask('>>> Cuts OK ? [y/n] ? [y] ')
                if not answ0:
                    answ0 = 'y'
                elif answ0 == 'no' or answ0 == 'NO':
                    answ0 = 'n'
        if goon:
            _z1, _z2 = str.split(str.split(sss[0])[0], '=')[1], str.split(str.split(sss[0])[1], '=')[1]
    else:
        print('Warning: image ' + str(img) + ' not found in the directory ')
    return _z1, _z2, goon


#############################################################################
def searchatmo(img, listatmo):

    hdr = floyds.util.readhdr(img)
    JD = readkey3(hdr, 'JD')
    _instrume = readkey3(hdr, 'TELID')
    camera = readkey3(hdr, 'INSTRUME')
    grism0 = readkey3(hdr, 'grism')
    if not listatmo:
        directory = floyds.__path__[0] + '/archive/' + str(_instrume) + '/' + camera +'/atmo/' + grism0
        listatmo = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listatmo:
        atmofile = ''
        distance = []
        goodlist = []
        for atmo in listatmo:
            print(atmo)
            hdra = readhdr(atmo)
            JDarc = readkey3(hdra, 'JD')
            grism1 = readkey3(hdra, 'grism')
            #           slit1=readkey3(hdra,'slit')
            if grism0 == grism1:
                goodlist.append(atmo)
                distance.append(np.abs(JD - JDarc))
        if len(distance) >= 1:
            atmofile = goodlist[np.argmin(distance)]
        else:
            atmofile = ''
    else:
        atmofile = ''
    return atmofile, directory


###########################################################################
def searcharc(img, listarc):

    hdr = floyds.util.readhdr(img)
    JD = readkey3(hdr, 'JD')
    _instrume = readkey3(hdr, 'TELID')
    camera = readkey3(hdr, 'INSTRUME')
    grism0 = readkey3(hdr, 'grism')
    slit0 = readkey3(hdr, 'slit')
    date0 = readkey3(hdr, 'date-night')
    #if slit0=='6.0' and _instrume in ['fts','2m0b']: slit0='2.0'
    #if slit0=='6.0' and _instrume in ['ftn','2m0a']: slit0='1.6'
    if not listarc:
        if str(_instrume) == 'fts' and camera == 'en12':
            in_era = datetime.datetime.strptime(str(date0), '%Y%m%d') > datetime.datetime(2024, 12, 1)
            in_era = in_era and datetime.datetime.strptime(str(date0), '%Y%m%d') < datetime.datetime(2025, 1, 15)
            if in_era:
                date_dir = 'post-2024-12-01/'
            else:
                date_dir = 'pre-2024-12-01/'
        else:
            date_dir = '' 
        directory = floyds.__path__[0] + '/archive/' + str(_instrume) + '/' + camera + '/arc/' + date_dir + grism0 + '/' + slit0
        listarc = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listarc:
        arcfile = ''
        distance = []
        goodlist = []
        for arc in listarc:
            print(arc)
            hdra = readhdr(arc)
            JDarc = readkey3(hdra, 'JD')
            grism1 = readkey3(hdra, 'grism')
            slit1 = readkey3(hdra, 'slit')
            if slit0 == slit1 and grism0 == grism1:
                goodlist.append(arc)
                distance.append(np.abs(JD - JDarc))
        if len(distance) >= 1:
            arcfile = goodlist[np.argmin(distance)]
        else:
            arcfile = ''
    else:
        arcfile = ''
    return arcfile, directory


###########################################################################
def searchsens(img, listsens):

    hdr = readhdr(img)
    if 'MJD' in hdr:
        JD = readkey3(hdr, 'MJD')
    elif 'MJD-OBS' in hdr:
        JD = readkey3(hdr, 'MJD-OBS')

    _instrume = readkey3(hdr, 'TELID')
    camera = readkey3(hdr, 'INSTRUME')
    grism0 = readkey3(hdr, 'grism')
    if not listsens:
        directory = floyds.__path__[0] + '/archive/' + str(_instrume) + '/' + camera + '/sens/' + grism0
        listsens = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listsens:
        sensfile = ''
        distance = []
        goodlist = []
        for sens in listsens:
            hdrs = readhdr(sens)
            if 'MJD' in hdrs:
                JDsens = readkey3(hdrs, 'MJD')
            elif 'MJD-OBS' in hdrs:
                JDsens = readkey3(hdrs, 'MJD-OBS')

            grism1 = readkey3(hdrs, 'grism') if readkey3(hdrs, 'grism') else readkey3(hdrs, 'GRISM')
            if grism0 == grism1:
                goodlist.append(sens)
                distance.append(np.abs(JD - JDsens))
        if len(distance) >= 1:
            sensfile = goodlist[np.argmin(distance)]
        else:
            sensfile = ''
    else:
        sensfile = ''
    return sensfile, directory


###########################################################################
def searchflat(img, listflat):

    hdr = readhdr(img)
    JD = readkey3(hdr, 'JD')
    _instrume = readkey3(hdr, 'TELID')
    camera = readkey3(hdr, 'INSTRUME')
    grism0 = readkey3(hdr, 'grism')
    if not listflat:
        directory = floyds.__path__[0] + '/archive/' + str(_instrume) + '/' + camera + '/flat/' + grism0
        listflat = glob.glob(directory + '/*fits')
    else:
        directory = ''
    if listflat:
        faltfile = ''
        distance = []
        goodlist = []
        for flat in listflat:
            hdrf = readhdr(flat)
            JDflat = readkey3(hdrf, 'JD')
            filter1 = readkey3(hdrf, 'filter')
            if filter0 == filter1:
                goodlist.append(flat)
                distance.append(np.abs(JD - JDflat))
        if len(distance) >= 1:
            flatfile = goodlist[np.argmin(distance)]
        else:
            flatfile = ''
    else:
        flatfile = ''
    return flatfile, directory


###########################################################################
def readstandard(standardfile):

    if os.path.isfile(standardfile):
        listastandard = standardfile
    elif standardfile[0] == '/':
        listastandard = standardfile
    else:
        listastandard = floyds.__path__[0] + '/standard/stdlist/' + standardfile
    f = open(listastandard, 'r')
    liststd = f.readlines()
    f.close()
    star, ra, dec = [], [], []
    magnitude = []
    for i in liststd:
        if i[0] != '#':
            star.append(str.split(i)[0])
            _ra = str.split(str.split(i)[1], ':')
            _dec = str.split(str.split(i)[2], ':')
            ra.append((float(_ra[0]) + ((float(_ra[1]) + (float(_ra[2]) / 60.)) / 60.)) * 15)
            if '-' in str(_dec[0]):
                dec.append((-1) * (np.abs(float(_dec[0])) + ((float(_dec[1]) + (float(_dec[2]) / 60.)) / 60.)))
            else:
                dec.append(float(_dec[0]) + ((float(_dec[1]) + (float(_dec[2]) / 60.)) / 60.))
            try:
                magnitude.append(str.split(i)[3])
            except:
                magnitude.append(999)
    return np.array(star), np.array(ra), np.array(dec), np.array(magnitude)

###########################################################################
def pval(_xx, p):
    _y = +p[0] + p[1] * _xx
    return _y


def residual(p, y, x):
    for i in range(len(p)):
        err = (y - p[i] * x ** i)
    return err

#########################################################################
def defsex(namefile):

    sexfile = floyds.__path__[0] + '/standard/sex/default.sex'
    f = open(sexfile, 'r')
    ss = f.readlines()
    f.close()
    ff = open(namefile, 'w')
    for i in ss:
        if string.count(i, 'PARAMETERS_NAME') == 1:
            ff.write('PARAMETERS_NAME  "' + floyds.__path__[0] + '/standard/sex/default.param"\n')
        elif string.count(i, 'FILTER_NAME') == 1:
            ff.write('FILTER_NAME  "' + floyds.__path__[0] + '/standard/sex/default.conv"\n')
        elif string.count(i, 'STARNNW_NAME') == 1:
            ff.write('STARNNW_NAME "' + floyds.__path__[0] + '/standard/sex/default.nnw"\n')
        else:
            ff.write(i)
    ff.close()
    return namefile


############################################################
def defswarp(namefile, imgname, _combine, gain=''):

    if _combine.lower() in ['median']:
        _combine = 'MEDIAN'
    elif _combine.lower() in ['average']:
        _combine = 'AVERAGE'
    elif _combine.lower() in ['sum']:
        _combine = 'SUM'
    swarpfile = floyds.__path__[0] + '/standard/sex/default.swarp'
    f = open(swarpfile, 'r')
    ss = f.readlines()
    f.close()
    ff = open(namefile, 'w')
    for i in ss:
        if string.count(i, 'IMAGEOUT_NAME') == 1:
            ff.write('IMAGEOUT_NAME    ' + str(imgname) + '  # Output filename \n')
        elif string.count(i, 'WEIGHTOUT_NAME') == 1:
            ff.write('WEIGHTOUT_NAME   ' + str(
                re.sub('.fits', '.weight.fits', imgname)) + '  # Output weight-map filename  \n')
        elif string.count(i, 'COMBINE_TYPE') == 1:
            ff.write('COMBINE_TYPE    ' + str(_combine) + '  # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CHI2 \n')
        elif string.count(i, 'GAIN_DEFAULT') == 1:
            if gain:
                ff.write('GAIN_DEFAULT    ' + str(gain) + '  # Default gain if no FITS keyword found \n')
            else:
                ff.write(i)
        else:
            ff.write(i)
    ff.close()
    return namefile


#################################################################################
def airmass(img, overwrite=True, _observatory='lasilla'):

    iraf.astutil(_doprint=0)
    hdr = readhdr(img)
    if readkey3(hdr, 'UTC'):
        _UT = (readkey3(hdr, 'UTC') + (readkey3(hdr, 'exptime') / 2)) / 3600
        _date = readkey3(hdr, 'date-obs')
        _date = _date[0:4] + '-' + _date[4:6] + '-' + _date[6:8]
        _RA = readkey3(hdr, 'RA') / 15
        _DEC = readkey3(hdr, 'DEC')
        f = file('airmass.txt', 'w')
        f.write('mst = mst ("' + str(_date) + '",' + str(_UT) + ', obsdb ("' + str(_observatory) + '", "longitude"))\n')
        f.write(
            'air = airmass (' + str(_RA) + ',' + str(_DEC) + ',mst, obsdb ("' + str(_observatory) + '", "latitude"))\n')
        f.write('print(air)\n')
        f.close()
        _air = iraf.astcalc(image=img, command="airmass.txt", Stdout=1)[0]
        try:
            _air = float(_air)
        except:
            _air = 999
        delete('airmass.txt')
        if overwrite and _air < 99.:
            floyds.util.updateheader(img, 0, {'AIRMASS': [_air, 'mean airmass computed with astcalc']})
    else:
        _air = ''
    return _air


####################################################################################
def dvex():
    dv = {}
    dv['line'] = {'red': 300, 'blu': 1000}
    dv['std'] = {'_t_order': 3, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance', '_nsum': 30, '_t_step': 10, '_t_nsum': 10, '_lower': -10, '_upper': 10,
                 '_b_sample': '-25:-15,15:25', '_resize': 'no', '_b_naver': -15}
    dv['obj'] = {'_t_order': 3, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance', '_nsum': 40, '_t_step': 10, '_t_nsum': 10, '_lower': -5, '_upper': 5,
                 '_b_sample': '-25:-15,15:25', '_resize': 'no', '_b_naver': -15}
    dv['agn'] = {'_t_order': 3, '_t_niter': 50, '_t_sample': '*', '_t_nlost': 20, '_width': 10, '_radius': 10,
                 '_weights': 'variance', '_nsum': 40, '_t_step': 10, '_t_nsum': 10, '_lower': -13, '_upper': 13,
                 '_b_sample': '-35:-20,20:35', '_resize': 'no', '_b_naver': -15}
#  -3, -15  
#  order of fit 
#  sigma clipping 
#  linear polynomial 
#  cleaning
#  sigma clipping 
    return dv


#################################################################################################################

def phase3header(img):
    from numpy import max, min, isfinite

    img_data = fits.open(img)[0].data
    hdr = readhdr(img)
    hedvec0 = {'DATAMIN': [float(np.min(img_data[np.isfinite(img_data)])), ''],
               'DATAMAX': [float(np.max(img_data[np.isfinite(img_data)])), ''],
               'ORIGIN': ['ESO', 'European Southern Observatory'],
               'PROCSOFT': ['floyds_' + str(floyds.__version__), 'pipeline version']}
    if readkey3(hdr, 'filter'):
        hedvec0['FILTER'] = [readkey3(hdr, 'filter'), 'Filter name']
    if readkey3(hdr, 'gain'):
        hedvec0['GAIN'] = [readkey3(hdr, 'gain'), 'Conversion from electrons to ADU']
    if readkey3(hdr, 'esoid'):
        hedvec0['OBID'] = [readkey3(hdr, 'esoid'), 'Observation block ID']
    if readkey3(hdr, 'esoprog'):
        hedvec0['PROG_ID'] = [readkey3(hdr, 'esoprog'), 'ESO program identification']
    if readkey3(hdr, 'tech'):
        hedvec0['OBSTECH'] = [readkey3(hdr, 'tech'), 'Observation technique']
    floyds.util.updateheader(img, 0, hedvec0)


#########################################################################################
def name_duplicate(img, nome, ext):  ###########################

    dimg = readkey3(readhdr(img), 'DATE-OBS')
    listafile = glob.glob(nome + '_?' + ext + '.fits') + glob.glob(nome + '_??' + ext + '.fits')
    if len(listafile) == 0:
        nome = nome + "_1" + ext + '.fits'
    else:
        date = []
        for l in listafile:
            date.append(readkey3(readhdr(l), 'DATE-OBS'))
        if dimg in date:
            nome = listafile[date.index(dimg)]
        else:
            n = 1
            while nome + '_' + str(n) + str(ext) + '.fits' in listafile:
                n = n + 1
            nome = nome + '_' + str(n) + str(ext) + '.fits'
    return nome


###############################################################################
def correctobject(img, coordinatefile):

    scal = np.pi / 180.
    std, rastd, decstd, magstd = readstandard(coordinatefile)
    img = re.sub('\n', '', img)
#    correctcard(img)
    hdr = readhdr(img)
    _ra = readkey3(hdr, 'RA')
    _dec = readkey3(hdr, 'DEC')
    dd = np.arccos(
        np.sin(_dec * scal) * np.sin(decstd * scal) + np.cos(_dec * scal) * np.cos(decstd * scal) * np.cos((_ra - rastd) * scal)) * (
             (180 / np.pi) * 3600)
    if min(dd) < 5200:
        floyds.util.updateheader(img, 0, {'OBJECT': [std[np.argmin(dd)], 'Original target.']})
        aa, bb, cc = rastd[np.argmin(dd)], decstd[np.argmin(dd)], std[np.argmin(dd)]
    else:
        aa, bb, cc = '', '', ''
    return aa, bb, cc


##################################################################################################
def archivingtar(outputlist, nametar):

    print('### making a tar with pre-reduced frames ........ please wait')
    stringa = ' '.join(outputlist)
    delete(nametar)
    os.system('tar -zcvf ' + nametar + ' ' + stringa)
    print('### tar file: ' + nametar + '\n')


#################################################################
def repstringinfile(filein, fileout, string1, string2):

    f = open(filein, 'r')
    ss = f.readlines()
    f.close()
    f = open(fileout, 'w')
    for n in range(len(ss)):
        if string1 in ss[n]:
            f.write(re.sub(string1, string2, ss[n]))
        else:
            f.write(ss[n])
    f.close()


###################################################

def StoN(img, ran=50):

    xx, yy = readspectrum(img)
    sntot = []
    xxmed = []
    for j in range(2, int((xx[-1] - xx[0]) / ran) - 2):
        aa, bb = xx[0] + j * ran - int(ran / 2.), xx[0] + j * ran + int(ran / 2.)
        ww = np.asarray([i for i in range(len(xx)) if ((xx[i] >= aa) & (xx[i] < bb) )])
        snr = np.average(yy[ww]) / np.sqrt((sum(((yy[ww] - np.average(yy[ww]))) ** 2)) / (len(ww) - 1))
        sntot.append(snr)
        xxmed.append(xx[0] + j * ran)
    return np.mean(sntot)


################################################
def spectraresolution(img):

    hdr = readhdr(img)
    _instrume = readkey3(hdr, 'instrume')
    _slit = readkey3(hdr, 'slit')
    _grism = readkey3(hdr, 'grism')
    risoluzioni = {}
    risoluzioni['floyds'] = {}
    risoluzioni['floyds']['blu', '6.0'] = 99
    risoluzioni['floyds']['blu', '2.0'] = 99
    risoluzioni['floyds']['blu', '1.6'] = 99
    risoluzioni['floyds']['red', '6.0'] = 99
    risoluzioni['floyds']['red', '2.0'] = 99
    risoluzioni['floyds']['red', '1.6'] = 99
    if _instrume in risoluzioni.keys():
        if (_grism, _slit) in risoluzioni[_instrume].keys():
            return risoluzioni[_instrume][_grism, _slit]
        else:
            return ''
    else:
        return ''


##################################################
def classifyfast(fitsfile, program='snid'):

    iraf.onedspec(_doprint=0)
    imgasci = re.sub('.fits', '.asci', fitsfile)
    floyds.util.delete(imgasci)
    iraf.onedspec.wspectext(fitsfile + '[*,1,1]', imgasci, header='no')

    if program == 'snid':
        print('\n######################\nclassify with snid\n')
        os.system('snid plot=0 iquery=0 inter=0 verbose=0 ' + imgasci)
        f = open(re.sub('.asci', '_snid.output', imgasci), 'r')
        ss = f.readlines()
        f.close()
        ss = ss[ss.index('### type fraction/redshift/age ###\n') + 2:ss.index(
            '### rlap-ordered template listings ###\n') - 2]
        bb = {}
        for i in ss:
            if string.split(i)[0] in ['Ia', 'Ia-norm', 'Ia-91T', 'Ia-91bg', 'Ia-csm', 'Ia-pec', 'Ib', 'Ib-norm',
                                      'Ib-pec', 'IIb', 'Ic', 'Ic-norm', 'Ic-pec', 'Ic-broad', 'II', 'II-pec', 'IIn',
                                      'IIP', 'IIL', 'NotSN', 'AGN', 'GAL', 'LBV']:
                bb[str.split(i)[0]] = {'all': str.split(i)[1:], 'frac': str.split(i)[2],
                                          'phase': str.split(i)[6], 'red': str.split(i)[3]}
        _type, _frac, _phase = [], [], []
        for ii in np.argsort(np.array([bb[i]['frac'] for i in bb.keys()], float))[::-1]:
            _type.append(bb.keys()[ii])
            _frac.append(float(bb[bb.keys()[ii]]['frac']))
            _phase.append(float(bb[bb.keys()[ii]]['phase']))
    elif program == 'superfit':
        print('classifiy with suprfit')
    elif program == 'gelato':
        print('gelato')
    else:
        print('warning: program not found')
    trigger = False
    if _type[0] not in ['AGN', 'NotSN', 'Gal']:
        if _phase[0] <= 0:
            trigger = True
        else:
            if 'pec' in _type[0]:
                trigger = True
            elif _type[0] in ['Ia-csm']:
                trigger = True
    print('\n#########################\n')
    print('\n   Type        %        phase   (most probable)')
    print('%7s\t%7s\t%7s' % (str(_type[0]), str(_frac[0]), str(_phase[0])))
    print('\n   Type        %        phase    (second most probable)')
    print('%7s\t%7s\t%7s' % (str(_type[1]), str(_frac[1]), str(_phase[1])))
    if trigger:
        print('\n##################\n INTERESTING SN !!!!\n ACTIVATE FOLLOW-UP WITH FULL LCOGT NETWORK !!!!\n\n')
    else:
        print('\n##################\n BORING SN ...... \n\n')
    #     print '\n We report that a spectrum of '+fitsfile+' was obtained robotically on Aug XX with the FLOYDS spectrograph '+
    #     ' at "Faulkes Telescope XXX". The spectrum (range 320-1000 nm) shows it to be a SN Ia roughly one week before maximum light, and is consistent with the host galaxy (CGCG 425-26) redshift of z=0.027. Classification was performed via supernova spectrum cross correlation using SNID (Blondin & Tonry, 2007, ApJ, 666, 1024).'
    return _type, _frac, _phase


##############################

################################################

def spectraresolution2(img0, ww=25):

    iraf.onedspec(_doprint=0)

    id = 'database/id' + re.sub('.fits', '', img0)
    img = re.sub('arc_', '', img0)
    data, hdr = fits.getdata(img0, 0, header=True)
    crvals = readkey3(hdr, 'CRVAL1')
    cds = readkey3(hdr, 'CD1_1')
    xx = np.arange(len(data))
    yy = data
    aa = crvals + (xx) * cds
    #   read identified lines from id file
    f = open(id, 'r')
    ss = f.readlines()
    f.close()
    indices = [i for i, x in enumerate(ss) if "begin" in x]
    if len(indices) <= 1:
        dd = ss[indices[0]:len(ss)]
    else:
        dd = ss[indices[0]:indices[1]]
    #         dd=ss[indices[-1]:len(ss)]
    start = [i for i, x in enumerate(dd) if "features" in x][0] + 1
    stop = [i for i, x in enumerate(dd) if "function" in x][0]
    ff = dd[start:stop]
    lines = []
    if len(ff) > 0:
        for i in ff:    lines.append(float(str.split(i)[2]))
        print(lines)
        lines = np.compress((aa[0] < np.array(lines)) & (np.array(lines) < aa[-1]), np.array(lines))
        cursor = ''
        yym = np.interp(lines - ww, aa, yy)
        yyp = np.interp(lines + ww, aa, yy)
        for i in range(0, len(lines)):
            cursor = cursor + str(lines[i] - ww) + '  ' + str(yym[i]) + '  1   k\n'
            cursor = cursor + str(lines[i] + ww) + '  ' + str(yyp[i]) + '  1   k\n'
        cursor = cursor + str(lines[i] + ww) + '  ' + str(yyp[i]) + '  1   q\n'
        ff = open('_cursor', 'w')
        ff.write(cursor)
        ff.close()

        aaa = iraf.noao.onedspec.bplot(img0, cursor='_cursor', spec2='', new_ima='', overwri='yes', Stdout=1)
        fw = []
        for i in aaa[1:]: fw.append(float(str.split(str.split(i, '=')[-1], 'k')[0]))
        floyds.util.delete('_cursor')
        res = (aa[0] + ((aa[-1] - aa[0]) / 2)) / np.mean(fw)
    else:
        res = 9999
    return res


##################################################

def spectraresolution3(img0, ww=25):

    iraf.onedspec(_doprint=0)

    img = re.sub('arc_', '', img0)
    data, hdr = fits.getdata(img0, 0, header=True)
    crvals = readkey3(hdr, 'CRVAL1')
    cds = readkey3(hdr, 'CD1_1')
    xx = np.arange(len(data))
    yy = data
    aa = crvals + (xx) * cds

    maxtab, mintab = floyds.util.peakdet(yy, np.median(yy) * 10)
    if len(maxtab) <= 0:      maxtab, mintab = floyds.util.peakdet(yy, np.median(yy))
    if len(maxtab) <= 0:
        lines = []
    else:
        lines0, b = zip(*maxtab)
        lines = crvals + (np.array(lines0)) * cds
    if len(lines) > 0:
        lines = np.compress((aa[0] < np.array(lines)) & (np.array(lines) < aa[-1]), np.array(lines))
        cursor = ''
        yym = np.interp(lines - ww, aa, yy)
        yyp = np.interp(lines + ww, aa, yy)
        for i in range(0, len(lines)):
            cursor = cursor + str(lines[i] - ww) + '  ' + str(yym[i]) + '  1   k\n'
            cursor = cursor + str(lines[i] + ww) + '  ' + str(yyp[i]) + '  1   k\n'
        cursor = cursor + str(lines[i] + ww) + '  ' + str(yyp[i]) + '  1   q\n'
        ff = open('_cursor', 'w')
        ff.write(cursor)
        ff.close()

        aaa = iraf.noao.onedspec.bplot(img0, cursor='_cursor', spec2='', new_ima='', overwri='yes', Stdout=1)
        fw = []
        for i in aaa[1:]: fw.append(float(str.split(str.split(i, '=')[-1], 'k')[0]))
        floyds.util.delete('_cursor')
        res = (aa[0] + ((aa[-1] - aa[0]) / 2)) / np.mean(fw)
    else:
        res = 9999
    return res


##################################################

def peakdet(v, delta, x=None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    Returns two arrays
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
    if x is None:
        x = np.arange(len(v))
    v = np.asarray(v)
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    mn, mx = np.inf, -np.inf
    mnpos, mxpos = np.nan, np.nan
    lookformax = True
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        if lookformax:
            if this < mx - delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
    return np.array(maxtab), np.array(mintab)

######################################################################

def mjdtoday():
  today = datetime.datetime.utcnow()
  mjd0  = datetime.datetime(1858,11,17)
  mjd = (today - mjd0).days
  return mjd

def to_safe_filename(unsafe_string):
    valid_filename_characters = "-_.(){letters}{numbers}".format(letters=string.ascii_letters, numbers=string.digits)
    return ''.join(character for character in unsafe_string if character in valid_filename_characters)
