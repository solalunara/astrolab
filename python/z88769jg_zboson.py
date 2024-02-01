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
from adjustText import adjust_text;

FILE_NAME = '../data.csv';
FILE_NAME_PREDICTED = "../data_predicted.csv";

PLANCK_CONST = 6.582119569e-25; #GeV*s

def luminosity( beta, velocity ):
    """calculates tully fisher relation for given rotational velocity

    Args:
        beta (list<float>): parameters
        velocity (float or NDArray): rotational velocity of the galaxy (adjusted for inclination)

    Returns:
        float or NDArray: luminosity of a galaxy with the given rotational velocity
    """
    return beta[ 1 ] * velocity**beta[ 0 ];

def linear( beta, velocity ):
    return beta[ 1 ] + velocity * beta[ 0 ];

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
        
    data = np.zeros( (0, 8) );
    names = [];
    SKIPPED_FIRST_LINE = False
    for line in input_file:
        if ( not SKIPPED_FIRST_LINE ):
            SKIPPED_FIRST_LINE = True;

        else:
            split_up = line.split( ',' );
            temp = np.array([float(split_up[1]), float(split_up[2]), float(split_up[3]), float(split_up[4]), float(split_up[5]), float(split_up[6]), float(split_up[7]), float(split_up[8])]);
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
    
    ax.errorbar( xdata, ydata, xerr=sigma_x, yerr=sigma_y, fmt='o', label="known distance galaxies" );
    ax.errorbar( pred_x, pred_y, yerr=pred_sigma_y / 2, fmt='o', label="predicted distance galaxies" );
    x_line = np.linspace( np.min( xdata ), np.max( xdata ), 1000000 );
    ax.plot( x_line, luminosity( result, x_line ), label="line of best fit" );
    
    texts = [ax.text(xdata[i], ydata[i], labels[i], ha='center', va='center', fontsize=8) for i in range(len(labels))]
    pred_texts = [ax.text(pred_x[i], pred_y[i], pred_labels[i], ha='center', va='center', fontsize=8) for i in range(len(pred_labels))]
    for i in range( len( pred_texts ) ):
        texts.append( pred_texts[ i ] );
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='red', lw=0.5))
    #adjust_text(pred_texts, arrowprops=dict(arrowstyle='->', color='red'))
    
    ax.set_title( 'Tully Fisher Relation' );
    ax.set_xlabel( 'Rotational Velocity km/s' );
    ax.set_ylabel( 'Luminosity Density W/Hz E+27' );
    ax.legend();
    plt.savefig( 'TullyFisher.png' );
    plt.show();
    
### --- MAIN CODE EXECUTION --- ###

## Read Data From Files ##
data = read_data( FILE_NAME );
data_predicted = read_data( FILE_NAME_PREDICTED );

xdata = data[ 0 ][ :, 0 ];
sigma_x = data[ 0 ][ :, 1 ];
ydata = data[ 0 ][ :, 2 ];
sigma_y = data[ 0 ][ :, 3 ];
dist = data[ 0 ][ :, 4 ];
dist_err = data[ 0 ][ :, 5 ];
vel = data[ 0 ][ :, 6 ];
vel_err = data[ 0 ][ :, 7 ];

pred_x = data_predicted[ 0 ][ :, 0 ];
pred_sigma_x = data_predicted[ 0 ][ :, 1 ];
pred_y = data_predicted[ 0 ][ :, 2 ];
pred_sigma_y = np.abs( data_predicted[ 0 ][ :, 3 ] );
pred_dist = data_predicted[ 0 ][ :, 4 ];
pred_dist_err = data_predicted[ 0 ][ :, 5 ];
pred_vel = data_predicted[ 0 ][ :, 6 ];
pred_vel_err = data_predicted[ 0 ][ :, 7 ];

odr_data = RealData( xdata, ydata, sigma_x, sigma_y  );
odr_model = Model( luminosity );
odr = ODR( odr_data, odr_model, beta0 = [ 6, 1e-15 ], maxit=100000 );
odr.set_job( fit_type=0 );
output = odr.run();

residuals = output.res_var;
parameters = output.beta;
sigmas = output.sd_beta / np.sqrt( output.res_var );


#res = sp.curve_fit( luminosity_curvefit, xdata, ydata, None, sigma_y, True );
#parameters = res[ 0 ];
#sigmas = np.sqrt( np.diag( res[ 1 ] ) );

## Plot Results ##
print( str( residuals ) );
print( str( parameters[ 0 ] ) + " +/- " + str( sigmas[ 0 ] ) )
print( str( parameters[ 1 ] ) + " +/- " + str( sigmas[ 1 ] ) )
plot_result( xdata, ydata, sigma_x, sigma_y, parameters, pred_x, pred_y, pred_sigma_x, pred_sigma_y, sigmas, data[ 1 ], data_predicted[ 1 ] );


fig = plt.figure();
ax = fig.add_subplot( 111 );

ax.errorbar( dist, vel, xerr=dist_err, yerr=vel_err, fmt='o', label="known distance galaxies" );
ax.errorbar( pred_dist, pred_vel, xerr=pred_dist_err / 2, yerr=pred_vel_err, fmt='o', label="unknown distance galaxies" );

odr_data = RealData( np.append( dist, pred_dist ), np.append( vel, pred_vel ), np.append( dist_err, pred_dist_err ), np.append( vel_err, pred_vel_err ) );
odr_model = Model( linear );
odr = ODR( odr_data, odr_model, beta0 = [ 4, 1 ] );
odr.set_job( fit_type=0 );
output = odr.run();

parameters = output.beta;
sigmas = output.sd_beta / np.sqrt( output.res_var );
residuals = output.res_var;

print( str( residuals ) );
print( str( parameters[ 0 ] ) + " +/- " + str( sigmas[ 0 ] ) )
print( str( parameters[ 1 ] ) + " +/- " + str( sigmas[ 1 ] ) )

x_line = np.linspace( np.min( np.append( dist, pred_dist ) ), np.max( dist ), 1000000 );
ax.plot( x_line, linear( parameters, x_line ), label="line of best fit" );

texts = [ax.text(dist[i], vel[i], data[ 1 ][i], ha='center', va='center', fontsize=8) for i in range(len(data[ 1 ]))]
pred_texts = [ax.text(pred_dist[i], pred_vel[i], data_predicted[ 1 ][i], ha='center', va='center', fontsize=8) for i in range(len(data_predicted[ 1 ]))]
for i in range( len( pred_texts ) ):
    texts.append( pred_texts[ i ] );
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='red', lw=0.5))


ax.set_title( 'Universal Expansion Relation' );
ax.set_xlabel( 'distance Mpc' );
ax.set_ylabel( 'velocity km/s' );
ax.legend();
plt.savefig( 'HubbleConst.png' );
plt.show();