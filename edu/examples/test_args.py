#!/usr/bin/env python
import re
import sys
import os
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description = """Parse the arguments and set to a dictionary containing 
            information about arguments, such
            as filelist, IsGitm, var (number), cut, diff (difference with
            other plots), movie, ext (for movie), rate (framerate for movie),
            tec, winds (plot with winds), alt (to plot), lat (to plot),
            lon (to plot), IsLog, and help (display help) """, usage= """gitm_plot_one_alt.py -var=N 
                -tec -winds -cut=alt,lat,lon 
                -alt=alt -lat=lat -lon=lon -alog 
                -help [*.bin or a file]""" , epilog= """ At end, enter -file and list the file you want to plot.
                This code should work with GITM files (*.bin) and Aether netCDF files (*nc)""")

    parser.add_argument('-var', metavar = 'var',  \
                        action='store_const', const = 0, default =15, \
                        help = '   -var=number : number is variable to plot')
    parser.add_argument('-cut', metavar = 'cut',  default ='alt',\
                        help = '   -cut=alt,lat,lon : which cut you would like')
    parser.add_argument('-wind', default =False,\
                        help='-winds: overplot winds', \
                        action="store_true")
    parser.add_argument('-alt', metavar = 'alt', default =400.0, \
                        help = '   -alt=altitude : can be either alt in km or grid number (closest)')
    parser.add_argument('-lat', metavar = 'lat',  default =-100.0,  \
                        help = '   -lat=latitude : latitude in degrees (closest)')
    parser.add_argument('-lon', metavar = 'lon',  default =-100.0,\
                        help = '   -lon=longitude: longitude in degrees (closest)')
    parser.add_argument('-alog',  default=False, action="store_true",\
                        help = '   -alog : plot the log of the variable')
    parser.add_argument('-diff', metavar = 'diff',  default =0, \
                        help = ' ')
    parser.add_argument('-mkv', metavar = 'mkv', action='store_const', const = "mkv",\
                        help = '   -mkv :')
    parser.add_argument('-mp4', metavar = 'mp4', action='store_const', const = 'mp4',\
                        help = '   -mp4 : ')
    parser.add_argument('-gif', metavar = 'gif', action='store_const', const = 'gif', \
                        help = '   -gif : ')
    parser.add_argument('-movie', metavar = 'movie',  default =0,\
                        action='store_const', const = '1', help = '   -movie : pe variable')
    parser.add_argument('-tec', metavar = 'tec',  default =0, \
                        help = '   -tec : ')
    parser.add_argument('-rate', metavar = 'rate',  default =30,\
                        help = '   -rate : ')
    parser.add_argument('-IsGitm', default=False, \
                        help='', \
                        action="store_true")
    parser.add_argument('-HasHeader', default =False,\
                        help='', \
                        action="store_true")
    parser.add_argument('-filename', nargs = 2, \
                        help = '   -file: for a new file not included')
    
    args = parser.parse_args()
    parser.print_help()
    return args


args = parse_args()


#if all above arguments are false or not read in, and sys.argv > 1, do stuff
#if IsFound==0 and not(arg==argv[0]):
 #               filelist.append(arg)
  #              m = re.match(r'(.*)bin',arg)
   #             if m:
    #                IsGitm = 1
     #               HasHeader = 0
      #              # check for a header file:
       #             checkFile = glob(m.group(1)+"header")
        #               if (len(checkFile[0]) > 1):
         #                   HasHeaders = 1
          #      else:
           #         IsGitm = 0