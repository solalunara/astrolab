#!/usr/bin/env python
""" assignment 2 zboson

This program fits a curve to the cross section of a Z0 Boson as a function
of its energy and calculates the uncertanties on the parameters.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
<http://www.gnu.org/licenses/>.
"""

__author__ = "Luna Greenberg"
__email__ = __contact__ = "jackgphysics@gmail.com"
__date__ = "2023/12/13"
__deprecated__ = False
__license__ = "GPLv3"
__status__ = "Production"
__version__ = "0.0.3"

import numpy as np;
import matplotlib.pyplot as plt;
import scipy.optimize as sp;
from scipy.odr import ODR, Model, Data, RealData;
import math;

FILE_NAME = '../data.csv';

PLANCK_CONST = 6.582119569e-25; #GeV*s

def luminosity( beta, velocity ):
    """calculates tully fisher relation for given rotational velocity

    Args:
        beta (list<float>): parameters
        velocity (float or NDArray): rotational velocity of the galaxy (adjusted for inclination)

    Returns:
        float or NDArray: luminosity of a galaxy with the given rotational velocity
    """
    if ( beta[ 0 ].ndim == 3 and velocity.ndim == 1 ):
        velocity = np.reshape( velocity, (len( velocity ),1,1) );
    return beta[ 0 ] + beta[ 1 ] * velocity;

def chi_squared( function, data, *params ):
    """
    Calculates the chi-squared value for an arbitrary function

    Args:
        function (lambda): predictive function f(x, ...) where ... are params
        data (array): 2d array, rows are x/y/sigma and cols are each data point
        params: various function parameters

    Returns:
        float: total chi-squared value
    """
    predicted = function( params, data[ 0 ] );
    observed = data[ 1 ];
    error = data[ 2 ];
    if ( predicted.ndim == 3 ):
        observed = np.reshape( observed, (len( observed ),1,1) );
        error = np.reshape( error, (len( error ),1,1) );
    return np.sum( ( ( predicted - observed ) / error )**2, 0 );
    
def read_data(file_name):
    """
    Reads in data given file name.

    Parameters
    ----------
    file_name : string

    Returns
    -------
    data : np.array of floats
        Data should be in format [x, y, uncertainty on y]
    """
    try:
        input_file = open( file_name, 'r' );
    except:
        print( 'Error locating ' + file_name + ', make sure your file is in the right place' );
        exit( 0 );
        
    data = np.zeros( (0, 4) );
    SKIPPED_FIRST_LINE = False
    for line in input_file:
        if ( not SKIPPED_FIRST_LINE ):
            SKIPPED_FIRST_LINE = True;

        else:
            split_up = line.split( ',' );
            try:
                temp = np.array([float(split_up[0]), float(split_up[1]), float(split_up[2]), float(split_up[3])]);
                for i in range( 0, 4 ):
                    if ( math.isnan( temp[ i ] ) or 0 == temp[ i ] ):
                        temp = None;
                if ( temp is None ):
                    continue;
            except:
                print( 'Bad data found! Make sure your data is numeric, >0, and has three comma-separated columns.' );
                continue;

            data = np.vstack( (data, temp) );

    input_file.close()

    return data

def plot_result( xdata, ydata, sigma_x, sigma_y, result, sigma_params ):
    """plot data and line of best fit

    Args:
        xdata (NDArray): x values for each point
        ydata (Array): y values for each point
        sigma_x (Array): errors on the x values
        sigma_y (Array): errors on the y values
        result (Array): result of curve fit, result[ 0 ] being the parameters of the fit
    """
    fig = plt.figure();

    ax = fig.add_subplot( 111 );

    ax.errorbar( xdata, ydata, xerr=sigma_x, yerr=sigma_y, fmt='o' );
    x_line = np.linspace( np.min( xdata ), np.max( xdata ), 1000000 );
    ax.plot( x_line, luminosity( result, x_line ) );
    
    ax.set_title( 'Tully Fisher Relation' );
    ax.set_xlabel( 'Rotational Velocity km/s' );
    ax.set_ylabel( 'Luminosity Density J' );
    plt.savefig( 'TullyFisher.png' );
    plt.show();
    
### --- MAIN CODE EXECUTION --- ###

## Read Data From Files ##
data = read_data( FILE_NAME );

xdata = np.log( data[ :, 0 ] );
sigma_x = data[ :, 1 ] / np.exp( xdata );
ydata = np.log( data[ :, 2 ] );
sigma_y = data[ :, 3 ] / np.exp( ydata );


odr_data = RealData( xdata, ydata, sigma_x, sigma_y  );
odr_model = Model( luminosity );
odr = ODR( odr_data, odr_model, beta0 = [ 1, 6 ] );
odr.set_job( fit_type=2 );
output = odr.run();

parameters = output.beta;
sigmas = output.sd_beta;

## Plot Results ##
print( str( parameters[ 0 ] ) + " +/- " + str( sigmas[ 0 ] ) )
print( str( parameters[ 1 ] ) + " +/- " + str( sigmas[ 1 ] ) )
plot_result( xdata, ydata, sigma_x, sigma_y, parameters, sigmas );