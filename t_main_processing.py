import json
import t_laads_tools
from config import *
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import numpy as np
from matplotlib import pyplot as plt
import datetime
import os

# Global dictionaries to make for efficient indexing and lookup
with open(f'{interrim_path}poly_info.json') as f:
    poly_info = json.load(f)

with open(f'{interrim_path}poly_to_tile_to_pixel.json') as f:
    poly_to_tile_to_pixel = json.load(f)

with open(f'{interrim_path}tile_to_poly.json') as f:
    tile_to_poly = json.load(f)

with open(f'{interrim_path}VNP46A2_laads_urls.json') as f:
    VNP46A2_laads_urls = json.load(f)


class VNP46A2Tile:

    def __init__(self, tile_id):

        self.tile_id = tile_id
        self.polygons = []


class FUAPolygon:

    def __init__(self, poly_id, tile_id):
        # ID of the polygon
        self.id = poly_id
        # Tile ID of the polygon (note: polygons may have more than one tile)
        self.tile_id = tile_id
        # Reference to previous Observation object
        self.prev_obs = None
        # Dictionary of results
        self.observations = []

    # Flatten and save the polygon data
    def save(self):
        # Output dictionary
        output_dict = {
                        'poly_id': self.id,
                        'tile_id': self.tile_id,
                        'observations': {}
                      }
        # For each observation
        for observation in self.observations:
            # Add year to dictionary if necessary
            if observation.year not in output_dict['observations'].keys():
                output_dict['observations'][observation.year] = {}
            # Add doy key and results
            output_dict['observations'][observation.year][observation.doy] = {}
            # For attribute in observation dictionary
            for attr in observation.__dict__.keys():
                # If the attribute is not the pixel dictionary, or the polygon reference
                if attr != 'px_dict' and attr != 'polygon':
                    # Transfer the attribute
                    output_dict['observations'][observation.year][observation.doy][attr] = getattr(observation, attr)
        # Open the output file
        with open(f'{output_path}{self.id}_{self.tile_id}.json', 'w', encoding="utf-8") as of:
            # Save the dictionary
            json.dump(output_dict, of, indent=4)


class Observation:

    def __init__(self, year, doy, dos, polygon):

        # Year of observation
        self.year = year
        # Day of Year of observation
        self.doy = doy
        # Day of Study of observation
        self.dos = dos
        # Reference to polygon object
        self.polygon = polygon
        # Total NTL
        self.total_ntl = 0
        # Total area of pixels
        self.total_area = 0
        # Total number of pixels
        self.total_pixels = 0
        # Quality flags
        # Filled pixels
        self.filled_pixels = 0
        # Poor quality pixels (QA Flag value 2 in Mandatory QA Flag)
        self.poor_quality_pixels = 0
        # Ephemeral pixels (QA Flag value 1 in Mandatory QA Flag)
        self.ephemeral_pixels = 0
        # Snow/Ice pixels (Flag value 1 in Snow/Ice Flag)
        self.snow_ice_pixels = 0
        # Membership metrics
        # Dictionary of pixels and their NTL (cleared after use)
        self.px_dict = {}
        # Count of pixels the same as previous observation
        self.px_same_as_prev = 0
        # Count of pixels lost since the previous observation
        self.px_lost_since_prev = 0
        # Count of pixels gained since the previous observation
        self.px_gained_since_prev = 0
        # Apples to Apples (total proportional change in NTL in pixels that were present in previous observation)
        self.px_same_as_prev_ntl_change = 0

    def compareToPrev(self):
        # Previous observation switch
        prev_obs = True
        # If there is a previous observation
        if self.polygon.prev_obs is not None:
            # If the previous observation is the previous day-of-study
            if self.polygon.prev_obs.dos == self.dos - 1:
                # For each pixel x in the previous observation
                for prev_x in self.polygon.prev_obs.px_dict.keys():
                    # For each pixel y in the previous observation
                    for prev_y in self.polygon.prev_obs.px_dict[prev_x].keys():
                        # Recover the previous NTL
                        prev_ntl = self.polygon.prev_obs.px_dict[prev_x][prev_y]
                        # If the pixel is in the current observation
                        if prev_x in self.px_dict.keys():
                            if prev_y in self.px_dict[prev_x].keys():
                                # Add pixels same as previous
                                self.px_same_as_prev += 1
                                # If the previous NTL is not 0
                                if prev_ntl != 0:
                                    # Add proportional change to the NTL change
                                    self.px_same_as_prev_ntl_change += (self.px_dict[prev_x][prev_y] - prev_ntl) / prev_ntl
                                # Otherwise (previous NTL is 0)
                                else:
                                    # Add the value (treat denominator as 1)
                                    self.px_same_as_prev_ntl_change += self.px_dict[prev_x][prev_y]
                # Calculate pixels gained since previous observation
                self.px_gained_since_prev = self.total_pixels - self.px_same_as_prev
                # Calculate pixels lost since previous observation
                self.px_lost_since_prev = self.polygon.prev_obs.total_pixels - self.px_same_as_prev
            # Otherwise
            else:
                # Switch to false
                prev_obs = False
                # Ditch the previous observation's dictionary
                self.polygon.prev_obs.px_dict = None
        # Otherwise
        else:
            # Switch to false
            prev_obs = False
        # If there was no viable previous observation
        if prev_obs is False:
            # Set same, gained, lost and ntl change since previous to None
            self.px_gained_since_prev = None
            self.px_lost_since_prev = None
            self.px_same_as_prev = None
            self.px_same_as_prev_ntl_change = None

# Process all the tiles (will resume if part-way through)
def process_all_tiles():
    # List of existing files
    file_list = []
    # Retrieve processed file details
    for root, dir, files in os.walk(output_path, topdown=False):
        # For each file name
        for name in files:
            # Append to list
            file_list.append(name)
    # List of tiles to process
    tile_list = []

    # For each tile
    for tile in tile_to_poly.keys():
        # For each polygon:
        for poly in tile_to_poly[tile]:
            # If the file doesn't exist
            if f'{poly}_{tile}.json' not in file_list:
                # If the tile is not already in the list
                if tile not in tile_list:
                    # Add tile to the list
                    tile_list.append(tile)
            # Otherwise
            else:
                # If the file is smaller than 1MB (stopped mid-writing or something)
                if os.stat(f'{output_path}{poly}_{tile}.json').st_size < 1E6:
                    # If the tile is not already in the list
                    if tile not in tile_list:
                        # Add tile to list for processing
                        tile_list.append(tile)
    # Send tile list to process
    process_multiple_tiles(tile_list)


def process_multiple_tiles(tile_list):

    # If the tile list is a single tile
    if isinstance(tile_list, list) is False:
        # Encapsulate in a list
        tile_list = [tile_list]

    # Create a LAADS session
    currSession = t_laads_tools.connect_to_laads()

    # Start a ProcessPoolExecutor
    process_exec = ProcessPoolExecutor(max_workers=5)
    # Submit the tiles
    for tile in tile_list:
        process_exec.submit(process_tile, tile, laads_session=currSession)
        print(f'Processing {tile}')


def process_tile(tile_id, laads_session=None):
    # If there is no existing session
    if laads_session is None:
        # Create a session
        laads_session = t_laads_tools.connect_to_laads()
    # Create tile object
    currTile = VNP46A2Tile(tile_id)
    # For each Polygon in the tile
    for poly_id in tile_to_poly[tile_id]:
        # Create a polygon object
        currTile.polygons.append(FUAPolygon(poly_id, tile_id))
    # For each year there are data for the tile
    for year in sorted(VNP46A2_laads_urls[tile_id].keys(), key=int):
        # For each doy there are data for the tile
        for doy in sorted(VNP46A2_laads_urls[tile_id][year].keys(), key=int):
            # Get the URL for the VNP46A2 file
            target_file_url = VNP46A2_laads_urls[tile_id][year][doy]
            # Assemble the full URL
            target_url = f"{VNP46A2_base_url}{year}/{doy}/{target_file_url}"
            # Get the H5 file on LAADS
            h5file = t_laads_tools.get_VNP46A2_file(laads_session, target_url)
            # While the status of the return is False (incomplete file?)
            while h5file is False:
                # Retry the request
                h5file = t_laads_tools.get_VNP46A2_file(laads_session, target_url)
            # Make numpy arrays of NTL data, Quality flags, and snow flags
            ntl_data = np.array(h5file['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['DNB_BRDF-Corrected_NTL'])
            qa_flags = np.array(h5file['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['Mandatory_Quality_Flag'])
            snow_flags = np.array(h5file['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['Snow_Flag'])
            # For each polygon
            for polygon in currTile.polygons:
                # Create an observation object
                currObs = Observation(int(year), int(doy), get_dos(year, doy), polygon)
                # Append the observation to the polygon
                polygon.observations.append(currObs)
                # For each pixel x
                for px_x in poly_to_tile_to_pixel[polygon.id][tile_id].keys():
                    # For each pixel y
                    for px_y in poly_to_tile_to_pixel[polygon.id][tile_id][px_x].keys():
                        # Get the NTL datum
                        ntl_datum = int(ntl_data[int(px_y)][int(px_x)])
                        # If the NTL datum is not a fill value
                        if ntl_datum != 65535:
                            # If the QA flag is not poor
                            if int(qa_flags[int(px_y)][int(px_x)]) != 2:
                                # If the pixel is not flagged for snow/ice:
                                if int(snow_flags[int(px_y)][int(px_x)]) != 1:
                                    # If the pixel is flagged for ephemerality (QA Flag value 1)
                                    if int(qa_flags[int(px_y)][int(px_x)]) == 1:
                                        # Add to count
                                        currObs.ephemeral_pixels += 1
                                    # Apply the 0.1 scale factor to convert to nW / cm^2 / sr
                                    ntl_datum *= 0.1
                                    # Add the pixel coordinates and value to pixel dictionary
                                    if px_x not in currObs.px_dict.keys():
                                        currObs.px_dict[px_x] = {}
                                    # Adjust the NTL value for the area of its pixel (weighting the sample by area)
                                    ntl_datum *= poly_to_tile_to_pixel[polygon.id][tile_id][px_x][px_y] * 1E10
                                    # Store the NTL value
                                    currObs.px_dict[px_x][px_y] = ntl_datum
                                    # Add to pixel count
                                    currObs.total_pixels += 1
                                    # Add to total NTL
                                    currObs.total_ntl += ntl_datum
                                    # Add to total area for the observation
                                    currObs.total_area += poly_to_tile_to_pixel[polygon.id][tile_id][px_x][px_y]
                                # Otherwise (snow/ice)
                                else:
                                    # Add to count
                                    currObs.snow_ice_pixels += 1
                            # Otherwise (bad QA Flag)
                            else:
                                # Add to poor QA flag count
                                currObs.poor_quality_pixels += 1
                                # Check if there was also a snow/ice flag
                                if int(snow_flags[int(px_y)][int(px_x)]) == 1:
                                    # Add to count
                                    currObs.snow_ice_pixels += 1
                            # Otherwise (filed value)
                        else:
                            currObs.filled_pixels += 1
                # Deal with the membership metrics etc. for the previous observation
                currObs.compareToPrev()
                # Set the polygon's previous observation to the current observation
                polygon.prev_obs = currObs
            # Print update
        print(f'Year {year} processed for tile {tile_id}')
    # For each polygon
    for polygon in currTile.polygons:
        # Save the dictionary for the polygon (polyid_tile.json)
        polygon.save()
    # Print update
    print(f'Finished processing {tile_id}')


# Return an ordered list of the tile files
def get_tile_file_list(tile_id, tile_dict):
    pass


# Get the day of study
def get_dos(year, doy):
    # Return day of study (since DOY 19, 2012)
    return ((int(year) - VNP46A2_first_year) * 365) + int(doy) + ((int(year) - VNP46A2_first_year) % 4) - VNP46A2_first_doy


# Optional visualization of the tile
def visualize_tile(tile_data, px_xs, px_ys):

    ntl_data = np.array(tile_data['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['DNB_BRDF-Corrected_NTL']).tolist()
    qa_flags = np.array(tile_data['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['Mandatory_Quality_Flag']).tolist()
    snow_flags = np.array(tile_data['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['Snow_Flag']).tolist()
    qa_flag_cloud = np.array(tile_data['HDFEOS']['GRIDS']['VNP_Grid_DNB']['Data Fields']['QF_Cloud_Mask']).tolist()

    qa_values = []
    snow_values = []
    cloud_values = []

    for ntl_row, qa_row, snow_row, cloud_row in zip(ntl_data, qa_flags, snow_flags, qa_flag_cloud):
        for ntl, qa, snow, cloud in zip(ntl_row, qa_row, snow_row, cloud_row):

            if snow == 1:
                if qa not in qa_values:
                    qa_values.append(qa)

            # if ntl != 65535:
            #     if qa not in qa_values:
            #         qa_values.append(qa)
            #     if snow not in snow_values:
            #         snow_values.append(snow)
            #     if cloud not in cloud_values:
            #         cloud_values.append(cloud)

    print(f'QA Values (when NTL is not Fill): {qa_values}')
    print(f'Snow Values (when NTL is not Fill): {snow_values}')
    print(f'Cloud Values (when NTL is not Fill): {cloud_values}')

    #plt.imshow(tile_data)
    #plt.scatter(px_xs, px_ys, 1, 'r', marker='s')

    plt.show()


# Get the average of NTL for a given window
def getWindowAverage(poly_dict, year, doy, window_size):
    # Get a list of day deltas
    indList = list(range(int(-np.ceil(window_size/2)), int(np.floor(window_size/2))))
    # Counter for observations used
    obsUsed = 0
    # Total NTL
    totalNTL = 0
    # For each day delta
    for delta in indList:
        # Get a datetime object for the target day
        currDT = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(doy) - 1) + datetime.timedelta(delta)
        targetDOY = currDT.timetuple().tm_yday
        # If the year is in the polygons observation years
        if str(currDT.year) in poly_dict['observations'].keys():
            # If the day is in the year subdict
            if str(targetDOY) in poly_dict['observations'][str(currDT.year)].keys():
                # Adjust the NTL by area
                ntl = float(poly_dict['observations'][str(currDT.year)][str(targetDOY)]['total_ntl'])
                area = float(poly_dict['observations'][str(currDT.year)][str(targetDOY)]['total_area']) * 1E10
                # If the NTL is not 0
                if ntl != 0 and area != 0:
                    # Add to NTL
                    totalNTL += ntl / area
                    # Add to the obs count
                    obsUsed += 1
    # If there was at least one observation used
    if obsUsed > 0:
        # Return the average
        return totalNTL / obsUsed
    # Otherwise (no observations)
    else:
        # Return 0
        return 0


# Get the month and the day from the year and the doy
def getMonthDay(year, doy):
    # Get a datetime object for the year/doy
    currDT = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(doy) - 1)
    # Return month and day of month
    return currDT.timetuple().tm_mon, currDT.timetuple().tm_mday


# Convert all files to ASCII
def convert_all_to_ascii():

    # For each polygon in the polygon > tile > pixel dictionary
    for polygon in poly_to_tile_to_pixel.keys():
        # If there is more than one tile (combined file expected)
        if len(poly_to_tile_to_pixel[polygon].keys()) > 1:
            # Load the json file
            with open(f'{output_path}{polygon}_combined.json', 'r') as f:
                # Load the polygon dictionary
                polygon_dict = json.load(f)
        # Otherwise (one tile file)
        else:
            # Get the tile name
            tile = list(poly_to_tile_to_pixel[polygon].keys())[0]
            # Load the json file
            with open(f'{output_path}{polygon}_{tile}.json', 'r') as f:
                # Load the polygon dictionary
                polygon_dict = json.load(f)
        # Get output filename
        output_file = f'{polygon}.txt'
        # Open file for writing
        f = open(f'{ascii_output_path}{output_file}', 'w')
        # Write the column titles
        f.write(
            'SITE MONTH DATE YEAR NTL AREA_(km^2) 7_DAY_AVERAGE 14_DAY_AVERAGE 30_DAY_AVERAGE TOTAL_PIXELS SAME_PIXELS GAINED_PIXELS LOST_PIXELS\n')
        # For each year (sorted)
        for year in sorted(polygon_dict['observations'].keys(), key=int):
            # For each day (sorted)
            for doy in sorted(polygon_dict['observations'][year].keys(), key=int):
                # Get the observation dictionary
                obs_dict = polygon_dict['observations'][year][doy]
                # Adjust the NTL by area
                ntl = float(obs_dict['total_ntl'])
                area = float(obs_dict['total_area'])

                if area != 0:
                    ascii_NTL = ntl / (area * 1E10)
                else:
                    ascii_NTL = 'NaN'

                ascii_area = area
                # Get the rolling averages
                ascii_7day = getWindowAverage(polygon_dict, year, doy, 7)

                if ascii_7day == 0:
                    ascii_7day = 'NaN'

                ascii_14day = getWindowAverage(polygon_dict, year, doy, 14)

                if ascii_14day == 0:
                    ascii_14day = 'NaN'

                ascii_30day = getWindowAverage(polygon_dict, year, doy, 30)

                if ascii_30day == 0:
                    ascii_30day = 'NaN'

                ascii_poly_id = polygon_dict['poly_id']
                ascii_month, ascii_dom = getMonthDay(year, doy)
                ascii_pixels = obs_dict['total_pixels']
                ascii_same = obs_dict['px_same_as_prev']

                if ascii_same is None:
                    ascii_same = 'NaN'

                ascii_gain = obs_dict['px_gained_since_prev']

                if ascii_gain is None:
                    ascii_gain = 'NaN'

                ascii_lost = obs_dict['px_lost_since_prev']

                if ascii_lost is None:
                    ascii_lost = 'NaN'

                # Write the data
                f.write(f'{ascii_poly_id} '
                        f'{ascii_month} '
                        f'{ascii_dom} '
                        f'{year} '
                        f'{ascii_NTL} '
                        f'{ascii_area} '
                        f'{ascii_7day} '
                        f'{ascii_14day} '
                        f'{ascii_30day} '
                        f'{ascii_pixels} '
                        f'{ascii_same} '
                        f'{ascii_gain} '
                        f'{ascii_lost}\n')
        # Close the file
        f.close()

# Converting the output from the process_tile function to an ASCII file
def convert_to_ascii(polygon_file):

    # Load the json file
    with open(f'{output_path}{polygon_file}', 'r') as f:
        polygon_dict = json.load(f)

    # Get output filename
    output_file = polygon_file.replace('.json', '.txt')

    # Open file for writing
    f = open(f'{ascii_output_path}{output_file}', 'w')
    # Write the column titles
    f.write('SITE MONTH DATE YEAR NTL AREA_(km^2) 7_DAY_AVERAGE 14_DAY_AVERAGE 30_DAY_AVERAGE TOTAL_PIXELS SAME_PIXELS GAINED_PIXELS LOST_PIXELS\n')
    # Reference the FUA polygon object
    #currPolygon = self.polygon_ids[polygon]
    # Get the dict of all years and days represented across the tiles
    #dateDict = currPolygon.getYearsDays()
    # For each year (sorted)
    for year in sorted(polygon_dict['observations'].keys(), key=int):
        # For each day (sorted)
        for doy in sorted(polygon_dict['observations'][year].keys(), key=int):
            # Get the observation dictionary
            obs_dict = polygon_dict['observations'][year][doy]
            # Adjust the NTL by area
            ntl = float(obs_dict['total_ntl'])
            area = float(obs_dict['total_area'])

            if area != 0:
                ascii_NTL = ntl / (area * 1E10)
            else:
                ascii_NTL = 0

            ascii_area = area
            # Get the rolling averages
            ascii_7day = getWindowAverage(polygon_dict, year, doy, 7)
            ascii_14day = getWindowAverage(polygon_dict, year, doy, 14)
            ascii_30day = getWindowAverage(polygon_dict, year, doy, 30)
            ascii_poly_id = polygon_dict['poly_id']
            ascii_month, ascii_dom = getMonthDay(year, doy)
            ascii_pixels = obs_dict['total_pixels']
            ascii_same = obs_dict['px_same_as_prev']
            ascii_gain = obs_dict['px_gained_since_prev']
            ascii_lost = obs_dict['px_lost_since_prev']
            # Write the data
            f.write(f'{ascii_poly_id} '
                    f'{ascii_month} '
                    f'{ascii_dom} '
                    f'{year} '
                    f'{ascii_NTL} '
                    f'{ascii_area} '
                    f'{ascii_7day} '
                    f'{ascii_14day} '
                    f'{ascii_30day} '
                    f'{ascii_pixels} '
                    f'{ascii_same} '
                    f'{ascii_gain} '
                    f'{ascii_lost}\n')
    # Close the file
    f.close()


# Combining per-tile files
def combine_per_tile_files():
    # For each polygon in the polygon > tile > pixel dictionary
    for polygon in poly_to_tile_to_pixel.keys():
        # If there is more than one tile
        if len(poly_to_tile_to_pixel[polygon].keys()) > 1:
            # Start a combined dictionary
            combined_dict = {}
            # Tile list
            tile_list = []
            # For each tile
            for tile in poly_to_tile_to_pixel[polygon].keys():
                # Add to tile list
                tile_list.append(tile)
                # Open the file
                with open(f'{output_path}{polygon}_{tile}.json', 'r') as f:
                    # Load the json
                    in_dict = json.load(f)
                # For each year in the observations
                for year in in_dict['observations'].keys():
                    # If the year is not in the combined dictionary yet
                    if year not in combined_dict.keys():
                        # Add the year
                        combined_dict[year] = {}
                    # For each day in the year
                    for day in in_dict['observations'][year].keys():
                        # If the day is not in the combined dictionary yet
                        if day not in combined_dict[year].keys():
                            # Add the day (and transfer the dictionary
                            combined_dict[year][day] = in_dict['observations'][year][day]
                        # Otherwise (day already there)
                        else:
                            # Add the total ntl
                            combined_dict[year][day]['total_ntl'] += in_dict['observations'][year][day]['total_ntl']
                            # Add the total area
                            combined_dict[year][day]['total_area'] += in_dict['observations'][year][day]['total_area']
                            # Add the total pixels
                            combined_dict[year][day]['total_pixels'] += in_dict['observations'][year][day]['total_pixels']
                            # Add the filled pixels
                            combined_dict[year][day]['filled_pixels'] += in_dict['observations'][year][day][
                                'filled_pixels']
                            # Add the poor quality pixels
                            combined_dict[year][day]['poor_quality_pixels'] += in_dict['observations'][year][day][
                                'poor_quality_pixels']
                            # Add the ephemeral pixels
                            combined_dict[year][day]['ephemeral_pixels'] += in_dict['observations'][year][day][
                                'ephemeral_pixels']
                            # Add the snow/ice pixels
                            combined_dict[year][day]['snow_ice_pixels'] += in_dict['observations'][year][day][
                                'snow_ice_pixels']
                            # If the incoming pixels same as previous is not None
                            if in_dict['observations'][year][day]['px_same_as_prev'] is not None:
                                # If the existing pixels same as previous is not None
                                if combined_dict[year][day]['px_same_as_prev'] is not None:
                                    # Add to all the totals
                                    combined_dict[year][day]['px_same_as_prev'] += in_dict['observations'][year][day]['px_same_as_prev']
                                    combined_dict[year][day]['px_lost_since_prev'] += in_dict['observations'][year][day]['px_lost_since_prev']
                                    combined_dict[year][day]['px_gained_since_prev'] += in_dict['observations'][year][day]['px_gained_since_prev']
                                    combined_dict[year][day]['px_same_as_prev_ntl_change'] += in_dict['observations'][year][day]['px_same_as_prev_ntl_change']
                                # Otherwise (existing is None, incoming is not None)
                                else:
                                    # Overwrite the Nones with the incoming
                                    combined_dict[year][day]['px_same_as_prev'] = in_dict['observations'][year][day]['px_same_as_prev']
                                    combined_dict[year][day]['px_lost_since_prev'] = in_dict['observations'][year][day]['px_lost_since_prev']
                                    combined_dict[year][day]['px_gained_since_prev'] = in_dict['observations'][year][day]['px_gained_since_prev']
                                    combined_dict[year][day]['px_same_as_prev_ntl_change'] = in_dict['observations'][year][day]['px_same_as_prev_ntl_change']
            # Encapsulate the combined dictionary with the general structure for the saved files
            outdict = { 'poly_id': in_dict['poly_id'],
                        'tile_id': tile_list,
                        'observations': combined_dict}
            # Save the combined file
            with open(f'{output_path}{polygon}_combined.json', 'w') as of:
                # Dump the dictionary
                json.dump(outdict, of, indent=4)