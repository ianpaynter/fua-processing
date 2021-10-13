from config import *
import requests
import json
from time import sleep
import h5py
import io

# Function to submit request to LAADS and keep trying until we get a response
def try_try_again(r, s, target_url):

    # If we get timed out
    while r.status_code != 200:
        # Print a warning
        print(f'Warning, bad response for {target_url}.')
        # Wait a hot second
        sleep(10)
        # Try again
        r = s.get(target_url)
    return r

# Connect to LAADS and return a session object
def connect_to_laads():
    # Header command utilizing security token
    authToken = {'Authorization': f'Bearer {laadsToken}'}
    # Create session
    s = requests.session()
    # Update header with authorization
    s.headers.update(authToken)
    # Return the session object
    return s

# Get a VNP46A2 H5 file from laads and return it as a numpy array
def get_VNP46A2_file(session_obj, target_url):

    # Request the H5 file from the provided URL
    r = session_obj.get(target_url)
    # If the request failed
    if r.status_code != 200:
        # Send to repeated submission function
        r = try_try_again(r, session_obj, target_url)
    # Try to convert into an h5 object
    try:
        # Convert to h5 file object
        h5file = h5py.File(io.BytesIO(r.content), 'r')
        # Convert the response content to an H5py File object and return
        return h5file
    # If it fails (incomplete file)
    except:
        # Print a warning
        print('Warning: File could not be converted to h5. Possibly incomplete.')
        # Return False
        return False

# Function to return a dictionary of VNP46A2 data available on LAADS
def get_VNP46A2_availability():

    # Data dictionary
    data_dict = {}
    # Header command utilizing security token
    authToken = {'Authorization': f'Bearer {laadsToken}'}
    # Target URL for laads data
    target_url = f"https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A2.json"
    # Create session
    s = requests.session()
    # Update header with authorization
    s.headers.update(authToken)
    # Get the years in json format from the target URL
    r = s.get(target_url)
    # Load the content of the response
    years = json.loads(r.text)
    # For each year in the data
    for year in years:
        # Get year value
        year_value = year["name"]
        # Construct year URL
        year_url = target_url.replace(".json", f"/{year_value}.json")
        # Get the days (adding the year to the original URL
        r = s.get(year_url)
        # If the request failed
        if r.status_code != 200:
            # Send to repeated submission function
            r = try_try_again(r, s, year_url)
        # Load the data as text
        days = json.loads(r.text)
        # For each day
        for day in days:
            # Get day value
            day_value = day["name"]
            # Construct day URL
            day_url = target_url.replace(".json", f"/{year_value}/{day_value}.json")
            # Get the tiles (adding the day and year to the URL)
            r = s.get(day_url)
            print(year_value, day_value)
            # If the request failed
            if r.status_code != 200:
                # Send to repeated submission function
                r = try_try_again(r, s, day_url)
            # Load the data as text
            tiles = json.loads(r.text)
            # For each of the tiles
            for tile in tiles:
                # Pull the file name of the tile
                file_name = tile["name"]
                # Split the name on the periods
                split_name = file_name.split('.')
                # Extract the year, day, and tile name
                tile_year = split_name[1][1:5]
                tile_doy = split_name[1][5:]
                tile_name = split_name[2]
                tile_collection = split_name[3]
                # Store the filename in the dictionary
                if tile_name not in data_dict.keys():
                    data_dict[tile_name] = {}
                if tile_year not in data_dict[tile_name].keys():
                    data_dict[tile_name][tile_year] = {}
                if tile_doy not in data_dict[tile_name][tile_year].keys():
                    data_dict[tile_name][tile_year][tile_doy] = file_name
    # Close the session
    s.close()

    return data_dict

# if __name__ == '__main__':
#
#     stime = time()
#
#     data_dict = get_VNP46A2_availability()
#
#     print(time() - stime)
#
#     dl_folder = Path(interrim_path)
#     filename = 'VNP46A2_laads_urls.json'
#     file = dl_folder.joinpath(filename)
#     # with open(file, 'wb') as of:
#     # of.write(r.content)
#     with open(file, 'w') as of:
#         json.dump(data_dict, of)




