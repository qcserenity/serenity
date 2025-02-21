# @file   showCube.py
#
# @date   Feb 19, 2025
# @author Anton Rikus
# @copyright \n
# This file is part of the program Serenity.\n\n
# Serenity is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.\n\n
# Serenity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.\n\n
# You should have received a copy of the GNU Lesser General
# Public License along with Serenity.
# If not, see <http://www.gnu.org/licenses/>.\n

import py3Dmol

def showCube(filename, iso=0.0025, element_colors=None, **kwargs):
    """
    Visualize a cube file using py3Dmol.

    Parameters:
    filename (str): Path to the cube file.
    iso (float): Isovalue for the volumetric data visualization. Default is 0.0025.
    element_colors (dict): Dictionary mapping element symbols to colors. Default is None.
    **kwargs: Additional keyword arguments for customization. Accepts 'vol_params', 'vol_params_negative', and 'styles'.

    Returns:
    None
    """
    # Load the cube file content
    with open(filename, 'r') as file:
        cube_data = file.read()

    v = py3Dmol.view()

    # Add volumetric data with customizable parameters
    vol_params = {
        'isoval': iso,
        'smoothness': 6,
        'opacity': 0.8,
        'voldata': cube_data,
        'volformat': 'cube.gz',
        'color': 'blue'
    }
    vol_params.update(kwargs.get('vol_params', {}))
    v.addVolumetricData(cube_data, "cube.gz", vol_params)

    vol_params_negative = vol_params.copy()
    vol_params_negative['isoval'] = -1.0 * iso
    vol_params_negative['color'] = 'red'
    vol_params_negative.update(kwargs.get('vol_params_negative', {}))
    v.addVolumetricData(cube_data, "cube.gz", vol_params_negative)

    v.addModel(cube_data, 'cube')

    # Define default colors for different elements
    default_element_colors = {
        'H': 'white',
        'C': 'grey',
        'F': 'pink',
    }
    if element_colors:
        default_element_colors.update(element_colors)

    # Apply styles based on element colors
    for element, color in default_element_colors.items():
        v.setStyle({'elem': element}, {'sphere': {'radius': 0.3, 'color': color}, 'stick': {'radius': 0.1, 'color': color}})

    # Apply additional styles if provided
    if 'styles' in kwargs:
        for style in kwargs['styles']:
            v.setStyle(style['spec'], style['style'])

    v.zoomTo()
    v.show()