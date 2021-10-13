These are the main modules for the processing of the Full Urban Area (FUA) polygons from the VIIRS DNB VNP46A2 data. You should be able to get a reasonable idea of how this was done from reading the comments, especially in the process_tile function in t_main_processing.py.

Some of the functions here might be generally useful. For example:

Surveying available VNP46A2 data in LAADS:

t_laads_tools.py
    get_VNP46A2_availability

Connecting to LAADS with a LAADS token (get your own):

t_laads_tools.py
    connect_to_laads

Get a specified VNP46A2 file (with some failure protection):

t_laads_tools.py
    get_VNP46A2_file

Using concurrent futures to hit LAADS with parallel requests for data (WARNING: Do not increase the number too high here or you will look like you're DDOS'ing LAADS):

t_main_processing.py
    process_multiple_tiles