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
FILE_NAME_PREDICTED = "../data_predicted.csv";

PLANCK_CONST = 6.582119569e-25; #GeV*s

def luminosity_curvefit( velocity, beta0, beta1 ):
    """calculates tully fisher relation for given rotational velocity

    Args:
        beta (list<float>): parameters
        velocity (float or NDArray): rotational velocity of the galaxy (adjusted for inclination)

    Returns:
        float or NDArray: luminosity of a galaxy with the given rotational velocity
    """
    if ( beta0.ndim == 3 and velocity.ndim == 1 ):
        velocity = np.reshape( velocity, (len( velocity ),1,1) );
    return np.log( beta0 ) + beta1 * velocity;

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
    predicted = function( data[ 0 ], *params );
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
    names = [];
    SKIPPED_FIRST_LINE = False
    for line in input_file:
        if ( not SKIPPED_FIRST_LINE ):
            SKIPPED_FIRST_LINE = True;

        else:
            split_up = line.split( ',' );
            temp = np.array([float(split_up[1]), float(split_up[2]), float(split_up[3]), float(split_up[4])]);
            names.append( split_up[ 0 ] );

            data = np.vstack( (data, temp) );

    input_file.close()

    return data, names

def plot_result( xdata, ydata, sigma_x, sigma_y, result, pred_x, pred_y, pred_sigma_x, pred_sigma_y, sigma_params, labels, pred_labels ):
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
    for i in range( 0, len( labels ) ):
        ax.annotate( labels[ i ], (xdata[ i ], ydata[ i ] + .1 * (-1)**(i) ) );
    ax.errorbar( pred_x, pred_y, xerr=pred_sigma_x, yerr=pred_sigma_y, fmt='o' );
    for i in range( 0, len( pred_labels ) ):
        ax.annotate( pred_labels[ i ], (pred_x[ i ], pred_y[ i ] + 1 * (-1)**(i) ) );
    x_line = np.linspace( np.min( xdata ), np.max( xdata ), 1000000 );
    ax.plot( x_line, luminosity_curvefit( x_line, *result ) );
    
    ax.set_title( 'Tully Fisher Relation' );
    ax.set_xlabel( 'Rotational Velocity km/s' );
    ax.set_ylabel( 'Luminosity Density J' );
    plt.savefig( 'TullyFisher.png' );
    plt.show();
    
### --- MAIN CODE EXECUTION --- ###

## Read Data From Files ##
data = read_data( FILE_NAME );
data_predicted = read_data( FILE_NAME_PREDICTED );

xdata = np.log( data[ 0 ][ :, 0 ] );
sigma_x = data[ 0 ][ :, 1 ] / np.exp( xdata );
ydata = np.log( data[ 0 ][ :, 2 ] );
sigma_y = data[ 0 ][ :, 3 ] / np.exp( ydata );

pred_x = np.log( data_predicted[ 0 ][ :, 0 ] );
pred_sigma_x = data_predicted[ 0 ][ :, 1 ] / np.exp( pred_x );
pred_y = np.log( data_predicted[ 0 ][ :, 2 ] );
pred_sigma_y = data_predicted[ 0 ][ :, 3 ] / np.exp( pred_y );

# odr_data = RealData( xdata, ydata, sigma_x, sigma_y  );
# odr_model = Model( luminosity );
# odr = ODR( odr_data, odr_model, beta0 = [ 1, 6 ] );
# odr.set_job( fit_type=2 );
# output = odr.run();

# parameters = output.beta;
# sigmas = output.sd_beta;

res = sp.curve_fit( luminosity_curvefit, xdata, ydata, None, sigma_y, True );
parameters = res[ 0 ];
sigmas = np.sqrt( np.diag( res[ 1 ] ) );

## Plot Results ##
print( str( parameters[ 0 ] ) + " +/- " + str( sigmas[ 0 ] ) )
print( str( parameters[ 1 ] ) + " +/- " + str( sigmas[ 1 ] ) )
plot_result( xdata, ydata, sigma_x, sigma_y, parameters, pred_x, pred_y, pred_sigma_x, pred_sigma_y, sigmas, data[ 1 ], data_predicted[ 1 ] );