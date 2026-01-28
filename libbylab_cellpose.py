print("\nLoading Packages ", end="...")
import argparse
from cellpose.io import imread, masks_flows_to_seg
from cellpose import models
from cellpose import version as cellpose_version
import cv2
from datetime import datetime
import numpy as np
import os
print(" Done!")

def find_edge_cells(masks, edge_buffer):
    """
    Find all cells that lie on the edge of the image and remove them.

    I could probably just check the bottom, top, left, and right edge for any
      ROIs that aren't 0 and remove all non-zero ROIs that lie on the 
      border
    """

    # Get list of ROIs within `edge_buffer` pixels of the edge of the image
    edge_ROIs = set(
        [0] + # Ensure the background (ROI=0) is excluded
        [   ROI 
            for list_of_ROIs in [
                set(masks[:,edge_buffer-1]),  # left edge
                set(masks[:,-1*edge_buffer]), # right edge
                set(masks[edge_buffer-1,:]),  # bottom edge
                set(masks[-1*edge_buffer,:]), # top edge
            ] 
            for ROI in list_of_ROIs 
        ]
    )
    
    return sorted(edge_ROIs)

def get_luminance_info(image, masks, ROIs):
    """
    Calculates the average brightness of cells and background
    
    Parameters
    ----------
    image : 2D array
        Raw image data. Each pixel value is of range 0-255
    masks: 2D array
        Cellpose output mask. Each pixel value corresponds to an ROI
    ROIs: list
        List containing ROI identifiers. ROI = Region of Interest
        ROI = 0 for background. ROI > 0 for cells
        
    Returns
    -------
    ROI_luminances: 2D array
        Array containing array of pixel luminances (pixel values) for each ROI
    """
    
    # Try to ensure data is usuable
    if isinstance(masks, list): 
        masks = np.array(masks)
    if len(masks.shape) != 2: 
        raise ValueError(f"Masks must be two dimensional but "
                         f"provided mask has {len(masks.shape)} dimensions.")
    
    ROI_luminances = []
    for ROI in ROIs:
        # Get coordinates of this ROI via the mask 
        coordinates = np.where(masks==ROI)
        x_coords = coordinates[1]
        y_coords = coordinates[0]

        # Get brightness levels for the identified coordinates of this cell
        luminances = []
        for x, y in zip(x_coords, y_coords):
            this_luminance = image[y][x]
            luminances.append(this_luminance)
        ROI_luminances.append(luminances)
    
    return ROI_luminances


def get_area_info(masks, ROIs):
    """
    Calculates the area of image consumed by cells and background
    
    Parameters
    ----------
    masks: 2D array
        Cellpose output mask. Each pixel value corresponds to an ROI
    ROIs: list
        List containing ROI identifiers. ROI = Region of Interest
        ROI = 0 for background. ROI > 0 for cells
        
    Returns
    -------
    area_of_image : int
        Area of the image in pixels
    area_of_ROIs : list
        List containing area (in pixels) for each ROI
    """
    
    # Try to ensure data is usuable
    if isinstance(masks, list): 
        masks = np.array(masks)
    if len(masks.shape) != 2: 
        raise ValueError(f"Masks must be two dimensional but "
                         f"provided mask has {len(masks.shape)} dimensions.")
    
    # Calculate Area
    area_of_image = masks.size
    area_of_ROIs = [ np.count_nonzero(masks==ROI) for ROI in ROIs ]
    
    return area_of_image, area_of_ROIs

def get_centers_widths_heights(masks, ROIs, plot_path, edge_buffer=1, save_images=False):
    """
    Calculates and returns the centers, widths, and heights of each ROI in the image.
    Each index in `centers`, `widths`, and `heights` corresponds to the ROI with that label.
    Any `centers[ROI]`, `widths[ROI]`, and `heights[ROI]` that are set to 0 were not sized
      (or had a size of 0).
    
    Parameters
    ----------
    masks: 2D array
        Cellpose output mask. Each pixel value corresponds to an ROI
    edge_buffer : int
        Number of pixels from the edge of the image an ROI must be to be sized
        
    Returns
    -------
    centers : numpy.ndarray
        Locations of the centers of the rectangles fit around each ROI in pixels.
    widths : list
        Width of each ROI in pixels
    heights : list
        Height of each ROI in pixels
    """

    # Initialize the arrays of widths and heights that will be returned
    centers = np.zeros((len(ROIs), 2))
    widths = np.zeros(len(ROIs))
    heights = np.zeros(len(ROIs))

    # Get width and height for each ROI and compile a photo of all boundaries
    total_image = np.zeros(masks.shape)
    for n_ROI, ROI in enumerate(ROIs):

        # Build an image with just this ROI
        ROI_locs = np.where(masks==ROI)
        ROI_image = np.zeros(masks.shape)
        ROI_image[ROI_locs] = 255
        ROI_image = ROI_image.astype(np.uint8)

        # Get the contouring of this ROI 
        contours, hierarchy = cv2.findContours(
            ROI_image, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        rectangle_center, (width, height), angle_of_rotation = cv2.minAreaRect(contours[0])

        # Save the center
        centers[n_ROI] = rectangle_center

        # Save the larger of width and height as the width so its always bigger
        if width >= height: 
            widths[n_ROI] = width
            heights[n_ROI] = height
        else: 
            widths[n_ROI] = height
            heights[n_ROI] = width

        # Construct the image of this ROI and its box then add it to the image
        rect = cv2.minAreaRect(contours[0])
        box = np.intp(cv2.boxPoints(rect))
        _ = cv2.drawContours(ROI_image, [box], 0, (36,255,12), 3)
        total_image += ROI_image

    # Plot all sized ROIs and their calculated boxes
    if save_images:
        fig, ax = matplotlib.pyplot.subplots()
        _ = ax.imshow(total_image)
        ax.set_title(f"After Removing ROIs within {edge_buffer} Pixels of the Edge")
        ax.set_xlabel("X [pixels]")
        ax.set_ylabel("Y [pixels]")
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig(plot_path)
        matplotlib.pyplot.close()
        del fig, ax

    return centers, widths, heights

def hist2d(
    x_array, y_array, 
    cmap="inferno", figsize=(5,4),
    bins=10, xyranges=None, 
    weights=None, density=None, 
    x_label=None, y_label=None, title=None,
    cbar_min=None, cbar_max=None, cbar_label=None, cbar_norm=None, cbar_scaling=1,
    save_name=None, white_background_limit=True
):
    """
    Plots a 2D histogram 
    """
    
    fig, ax = matplotlib.pyplot.subplots(figsize=figsize)
    if isinstance(cmap, matplotlib.colors.ListedColormap) and white_background_limit:
        cmap.set_under("w", white_background_limit)
    H, xbins, ybins = np.histogram2d(
        x_array, y_array, bins=bins, range=xyranges,
        weights=weights, density=density
    )
    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
    mappable = ax.imshow(H.T*cbar_scaling, extent=extent, cmap=cmap, 
                         origin='upper', aspect='equal', norm=cbar_norm)
    
    if x_label  !=None: ax.set_xlabel(x_label)
    if y_label  !=None: ax.set_ylabel(y_label)
    if title    !=None: ax.set_title(title)
    if cbar_min !=None: mappable.norm.vmin = cbar_min
    if cbar_max !=None: mappable.norm.vmax = cbar_max
    
    if cbar_label is None: 
        if isinstance(weights, np.ndarray) or isinstance(weights, list):
            cbar_label="Weighted Number of Events"
        else:
            cbar_label="Number of Events"
    matplotlib.pyplot.colorbar(mappable=mappable, label=cbar_label)

    matplotlib.pyplot.tight_layout()
    
    if save_name!=None: matplotlib.pyplot.savefig(save_name, dpi=300)
        
    return fig, ax

def compare_csv_files(old_info, updated_info):
    """
    about 

    Parameters
    ----------
    old_info : 2d array
        Same info as update_info but without first column
    updated_info : 2d array
        Contains cell segmenting infomaiton with the columns defined below
    """

    columns = [ "Location of Image", "Number of Cells", 
                "Mean Bkg Luminance", "Mean Cell Luminance",
                "Background Area", "Total Cell Area", "Mean Cell Area" ]
    for i in range(1, len(columns)):
        old_data     = old_info[:,i-1]
        updated_data = [ arr[i] for arr in updated_info ]

        old_avg        = int( np.mean(old_data) )
        old_stddev     = int( np.std(old_data) )
        updated_avg    = int( np.mean(updated_data) )
        updated_stddev = int( np.std(updated_data) )

        if i==1: 
            print(f"{columns[i]}: \t{old_avg} -->> {updated_avg}")
        else:
            print(f"{columns[i]}: \t{old_avg} +/- {old_stddev} -->> "
                f"{updated_avg} +/- {updated_stddev}")
    print()

    return


def main():
    """
    This function contains the meat of this analysis script. 
    This function oversees the handling of user input, segmenting of cells, 
      and saving relevant information into a csv file
    """
    
    # Read command line arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'dir_to_analyze', type=str,
        help="Path to folder containing images to analyze.")
    parser.add_argument(
        'model_type', type=str, default="cyto",
        help="Model to use. Give Cellpose model or path to your own.")
    parser.add_argument(
        '--cytoplasm_color', type=int, default=0,
        help="Color of cytoplasm. 0=Gray, 1=Red, 2=Green, 3=Blue")
    parser.add_argument(
        '--nucleus_color', type=int, default=0,
        help="Color of nucleus. 0=Gray, 1=Red, 2=Green, 3=Blue")
    parser.add_argument(
        '--cell_diameter', type=int, default=60,
        help="Expected diameter of cells, in pixels.")
    parser.add_argument(
        '--um_per_pixel', type=float, default=1.,
        help="Micrometers per pixel.")
    parser.add_argument(
        '--update_csv', type=str, default=None,
        help="If provided, code will update csv output in --dir_to_analyze.")
    parser.add_argument(
        '--edge_buffer', type=int, default=5,
        help=("ROIs within this many pixels of the edge are ignored in the "
              "area, width, and height calculation.")
    )
    parser.add_argument(
        '--image_extensions', type=str, nargs="*", default=["JPG"],
        help="Extensions of the images you intend to analyze. Case sensitive.")
    parser.add_argument(
        '--save_boundary_images', action="store_true",
        help=("If present, will save pngs depicting rectangular boundaries "
              "used to calculate cell area, width, and height.")
    )
    parser.add_argument(
        '--plot_isodensity_map', action="store_true",
        help=("If present, will save an isodensity map of the identified cells.")
    )
    parser.add_argument(
        '--isodensity_map_binsize', type=float, default=1,
        help=("Size of isodensity map bins in millimeters.")
    )
    args = parser.parse_args()
    
    # Get list of files in folder
    if args.dir_to_analyze[-1] != "/":
        args.dir_to_analyze += "/"
    dir_contents = os.listdir(args.dir_to_analyze)
    image_files = [ args.dir_to_analyze+f 
                    for f in dir_contents 
                    if f.split(".")[-1] in args.image_extensions]
    print(f"\nStarting script. Will analyze {len(image_files)} "
            f"images from:\n\t{args.dir_to_analyze}")
    image_files.sort()

    # Load plotting package if saving cell boundary images or plotting the isodensity maps
    if args.save_boundary_images or args.plot_isodensity_map: 
        globals()["matplotlib"] = __import__(
            "matplotlib", globals(), locals(), ['pyplot'], 0)
    
    # Setup Cellpose Model (if not updating)
    if args.update_csv == None:
        print(f"Using the model:\n\t{args.model_type}\n")
        if int(cellpose_version.split(".")[0]) >= 4: 
            model = models.CellposeModel(pretrained_model=args.model_type)
        else:
            if (
                not os.path.exists( os.path.join(str(models._MODEL_DIR_ENV), args.model_type) )
                and
                not os.path.exists( os.path.join(str(models._MODEL_DIR_DEFAULT), args.model_type) )
            ):
                print(
                    "WARNING: The model you passed is not available in one of",
                    "the places Cellpose looks for models in:\n        ",
                    [str(path) for path in [models._MODEL_DIR_ENV, models._MODEL_DIR_DEFAULT] if path is not None],
                    "\n         Cellpose may load a default model or crash while trying",
                    "to load the model you requested:\n        ", 
                    args.model_type,
                    "\n         You may need to move your model to one of the",
                    "expected locations listed above and try again.\n"
                )
            model = models.CellposeModel(model_type=args.model_type)
        if args.model_type != os.path.basename(str(model.pretrained_model[0])): 
            print(
                f"The model you requested ({args.model_type}) is different from "
                f"the one Cellpose loaded ({model.pretrained_model})"
            )
            if input("Would you like to continue? (y/n)").strip().lower() != "y":
                print("Ending script")
                return -1
        channels = [args.cytoplasm_color, args.nucleus_color]
    
    # Load old csv file and updated seg files (if we're updating csv)
    if args.update_csv != None:
        # Load csv to update and its comments
        with open(args.update_csv) as old_csv:
            comments = [next(old_csv) for _ in range(2)]
        model = comments[0][13:-1]
        channels = [int(x.strip()) for x in comments[1][12:].split("Cell diameter")[0].split(" ")]
        cell_diameter = float(comments[1].split(":")[-1][:-1])
        old_image_files = np.loadtxt(
            args.update_csv, skiprows=3, dtype=str, delimiter=",", usecols=[0])
        old_info_to_save = np.loadtxt(
            args.update_csv, skiprows=3, delimiter=",", usecols=[1,2,3,4,5,6]
        ).astype(np.int64)

        # Find all the seg files 
        seg_files = [ args.dir_to_analyze+f 
                      for f in dir_contents if f[-4:]==".npy" ]
        seg_files.sort()
        print(f"\nStarting script. Will update cellpose_output.csv for "
              f"{len(seg_files)} images from:\n\t{args.dir_to_analyze}")
    
    # Run segmenter
    update_every = 24
    info_to_save = []   # Collects information about to be saved to csv file
    start_time = datetime.now()
    for i in range(len(image_files)):
        # Tell user how many files have been analyzed
        if i%update_every == 0 or i<5: 
            print(f"{i} files analyzed", end="")
        else: 
            print(i, end="", flush=True)
        
        # Get image and load with Cellpose
        cellpose_image = imread(image_files[i])
        print(".", end="", flush=True)
        
        # Segment the cells and save Cellpose output to seg file
        if args.update_csv != None:
            seg_file = np.load(seg_files[i], allow_pickle=True)[()]
            masks = seg_file['masks']
        else:
            masks, flows, styles = model.eval(
                cellpose_image, diameter=args.cell_diameter)
            print(".", end="", flush=True)
            if int(cellpose_version.split(".")[0]) >= 4: 
                masks_flows_to_seg(cellpose_image, masks, flows, image_files[i])
            elif int(cellpose_version.split(".")[0]) == 3: 
                masks_flows_to_seg(cellpose_image, masks, flows, image_files[i], 
                                   diams=args.cell_diameter, channels=channels)
            elif int(cellpose_version.split(".")[0]) == 2: 
                masks_flows_to_seg(cellpose_image, masks, flows, args.cell_diameter, 
                                   image_files[i], channels)
            else:
                raise NotImplementedError(
                    f"This script does not support Cellpose version {cellpose_version}")
            print(".", end="", flush=True)
        
        # Analyze the results
        ROIs = np.arange( 0, max(masks.flatten())+1 )
        edge_ROIs = find_edge_cells(masks, args.edge_buffer)
        ROIs_not_edge = [ROI for ROI in ROIs if ROI not in edge_ROIs]

        area_of_image, area_of_ROIs = get_area_info(masks, ROIs)
        area_of_ROIs_not_edge = np.array(area_of_ROIs)[ROIs_not_edge]
        luminances_of_ROIs_by_ROI = get_luminance_info(cellpose_image, 
                                                       masks, ROIs)
        luminances_of_all_ROIs = [ 
            L for ROI_Ls in luminances_of_ROIs_by_ROI 
            for L in ROI_Ls 
        ] # Flatten luminances array so it's compatible with numpy
        centers, widths, heights = get_centers_widths_heights(
            masks, ROIs_not_edge,
            image_files[i][ : image_files[i].rindex(".") ]+".pdf", 
            args.edge_buffer, 
            save_images=args.save_boundary_images
        )

        # Create the isodensity map, if one was requested
        # Make figure size real

        if args.plot_isodensity_map:
            millimeters_per_pixel = args.um_per_pixel/1000
            xs = centers[:,0] * millimeters_per_pixel # Convert location of cell centers to locations in millimeters
            ys = centers[:,1] * millimeters_per_pixel # Convert location of cell centers to locations in millimeters
            image_width = cellpose_image.shape[1] * millimeters_per_pixel # in millimeters
            image_height = cellpose_image.shape[0] * millimeters_per_pixel # in millimeters
            binsize = args.isodensity_map_binsize
            hist2d(
                xs, ys, 
                bins=(
                    np.arange( 0, ((image_width //binsize) + 1) * binsize, binsize), # X bins
                    np.arange( 0, ((image_height//binsize) + 1) * binsize, binsize) # Y bins
                ), # Get the number of bins in each dimension based on the image size and requested binsize
                figsize=(
                    4.8 * (image_width/image_height) + 1, # Width of final image (in inches) 
                        # size of image +1 inch for colorbar
                    4.8 # height of final image (in inches)
                ),
                x_label=r"X [mm]", y_label=r"Y [mm]", 
                cbar_label=r"Cells per mm$^2$", cbar_scaling=1/args.isodensity_map_binsize**2,
                save_name = image_files[i][ : image_files[i].rindex(".") ]+"_isodensity.png"
            )

        # Collect information to save
        number_of_cells = len(ROIs) - 1   # Background is ROI 0, subtract 1
        mean_background_luminance = np.mean(luminances_of_ROIs_by_ROI[0])
        mean_cell_luminance = np.mean(luminances_of_all_ROIs)
        background_area = area_of_ROIs[0] * (args.um_per_pixel**2)
        total_cell_area = np.sum(area_of_ROIs_not_edge[1:]) * (args.um_per_pixel**2)
        mean_cell_area = np.mean(area_of_ROIs_not_edge[1:]) * (args.um_per_pixel**2)
        mean_cell_width = np.mean(widths) * args.um_per_pixel
        mean_cell_height = np.mean(heights) * args.um_per_pixel
        info_to_save.append([
            image_files[i], 
            number_of_cells,
            mean_background_luminance, mean_cell_luminance, 
            background_area, total_cell_area, mean_cell_area,
            mean_cell_width, mean_cell_height
        ])
        
        # Tell user quantitative info about every [update_every] image
        # print (datetime.now()-start_time).seconds/60
        time_elapsed = (datetime.now()-start_time).seconds
        time_text = f"{time_elapsed:.0f} seconds" if time_elapsed<60 \
                    else f"{time_elapsed/60:.1f} minutes" 
        if i%update_every == 0 or i<5 : 
            print(
                f"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
                f"About image {image_files[i].split('/')[-1]}\n"
                f"\t{number_of_cells} cells identified.\n"
                f"\tAverage Luminance:"
                f"    Background: {mean_background_luminance:.0f}"
                f"    Cells: {mean_cell_luminance:.0f}\n"
                f"\tArea:    Background: {background_area:.0f}\n"
                f"\t\t Cells (total): {total_cell_area:.0f}"
                f"    Cells (average): {np.mean(mean_cell_area):.0f}\n"
                f"\tAverage Cell Width: {mean_cell_width:.0f}"
                f"    Average Cell Height: {mean_cell_height:.0f}\n"
                f"Total time elapsed: {time_text}\n"
                f"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            )
            
    # Tell user that we finished with the modeling 
    print(f"\n\nScript finished. "
          f"Total time elapsed: {(datetime.now()-start_time).seconds/60:.1f} "
          f"minutes.\n")

    # # If we're updating the csv, compare the old one to the new one
    # if args.update_csv != None:
    #     print("Comparing old csv to new csv.")
    #     compare_csv_files(old_info_to_save, info_to_save)
    
    # Save the info_to_save object to csv file
    if args.update_csv != None:
        date = datetime.now().strftime("%m%d%y")
        csv_file_path = args.dir_to_analyze+f"cellpose_output-{date}.csv"
        print(f"Saving new cellpose output to:\n{csv_file_path}\n")
    else:
        csv_file_path = args.dir_to_analyze+"cellpose_output.csv"
    header = (
        "Image Path, Number of cells, "
        "Average Luminance of Background, Average Luminance of Cells, "
        "Area of Background, Area of All Cells, Average Area of Single Cells,"
        "Average Cell Width, Average Cell Height"
    )
    comments = (
        f"# Model Used: {args.model_type}\n"
        f"# Channels: {' '.join([str(c) for c in channels])}\tCell diameter: {args.cell_diameter}\n"
    )
    np.savetxt(csv_file_path, info_to_save, delimiter=",", header=header, 
               fmt='% s', comments=comments)
    
    return 0

            
if __name__=="__main__": 
    main()