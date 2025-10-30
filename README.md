# LibbyLabCellpose

This repository contains the script that uses the [Cellpose](https://cellpose.readthedocs.io/en/latest/) module and models trained to identify retinal ganglion cells for analysis, first made for use in Professor [Richard Libby's lab](https://www.urmc.rochester.edu/eye-institute/research/labs/libby) at the University of Rochester Medical Center. This synthesizes Cellpose outputs for one or more images into one summary comma-separated value file, compiling both new calculations and pre-existing analysis metrics for each file. This work was led by Sarah Yablonski (University of Rochester) and Abigail Bishop (University of Wisconsin Madison). Please contact Professor Richard Libby with questions and comments.

# Contents

- `libbylab_cellpose.py` is the script that takes a directory full of images, analyzes them, and returns a comma-separated value file with each image's average background luminance, average cell luminance, background area, total cell area, average cell area, average cell width, and average cell height.
- `models/` contains a variety of models for use with Cellpose curated to identify retinal ganglion cells (RGCs). 
    - `Brn3a_40X_final`: Annotates RGCs identified using the marker [BRN3A](https://pubmed.ncbi.nlm.nih.gov/19264888/), with images taken at 40X magnification
    - `RBPMS_10X_Keyence_final`: Annotates RGCs identified using [RBPMS](https://pubmed.ncbi.nlm.nih.gov/24318667/) in a full retina. This model was trained using images taken on a Keyence BZ-X800 epifluorescent microscope at 10X magnification.
    - `RBPMS_40X_final`: Annotates RGCs identified using the marker [RBPMS](https://pubmed.ncbi.nlm.nih.gov/24318667/), with images taken at 40X magnification
    - `TUJ1_40X_final`: Annotates RCGs identified using the marker [TUJ1](https://onlinelibrary.wiley.com/doi/abs/10.1002/neu.480220109), with images taken at 40X magnification

# Summary of Cellpose procedure:

1. Build your model in the Cellpose GUI
1. Run the `libbylab_cellpose.py` over a directory containing images to be analyzed
1. Open the images in the Cellpose GUI and make corrections
1. Run the script over the new npy files to update the csv summary file

# How to run the script initially:

You can do a lot of this by dragging and dropping files into your terminal. You shouldn't need to enter the directory containing the script with your Terminal if you use drag-and-drop. Drag-and-dropping should automatically add spaces between each item you drag into your Terminal but please be mindful if it doesn't.  

1. Open your terminal. If you have an environment, make sure it is initialized.
1. Type `python`
1. Type the path or drag-and-drop the libbylab_cellpose.py python file into your Terminal
1. Type the path or drag-and-drop the folder (containing all the images you want to analyze) into your Terminal
1. Type the path or drag-and-drop the model you want to use into your Terminal
1. If you want to change the expected cytoplasm color, nucleus color (both default to gray) or cell diameter (defaults to 60 pixels) you can specify this now. For example, type `--cytoplasm_color 2` for green cytoplasm (note the space before 2)

# What the script outputs: 

* `cellpose_output.csv` 
    * A summary of the analysis of each sample including the number of identified cells, the average luminance of the background and cells, the area of the background, the combined area of all cells, the average area of individual cells, the average cell width, and the average cell height. 
* `*.npy` files
    * Python output files that contain information about the cell identification in the sample. 
* `*png` files
	* An image of the sample that was analyzed with boxes around each cell depicting how the height and width was calculated. Only output if you run the script with the command line argument `-save_boundary_images`


# To run the script over new npy files (update the csv)

1. Open Terminal. If you have an environment, make sure it is initialized.
1. Type `python`
1. Type the path or drag-and-drop the libbylab_cellpose.py python file into your Terminal
1. Type the path or drag-and-drop the folder (containing all the images you analyzed and their new seg.npy files) into your Terminal
1. Type the path or drag-and-drop the model you used into your Terminal
1. Type `--update_csv ` (don't forget the space)
1. Type the path or drag-and-drop the csv file you'd like to update into your Terminal

# Standard Output while Running the Script

There is a lot of information printed to the screen so that you can monitor the script as it runs. If something seems amiss, you can quit the script by clicking the `Control`/`ctrl` key and the `C` key at the same time. Please make sure Terminal is your active window when doing this. 

The start of the script will tell you how many (JPG) images are in the folder you provided. It looks like this: 
```
Starting script. Will analyze 144 images from:
    /Users/abigailbishop/Documents/cellpose/
Using the model:
    /Users/abigailbishop/Documents/cellpose/RGC_Count
```

The script will tell you which image its analyzing. It only outputs the index of the image, but for the first 5 images and for every 24th image the script will say "x images analyzed"

Note: This output begins counting at 0, not at 1. So if you see 13 output in this context, you are actually in the process of analyzing your 14th image, not the 13th. That's why you see `About image ...` information for images 0-4 and not for images 1-5.

For each image, the script outputs a period 3 times (mimicking a progress bar)
1. After an image is loaded
1. After Cellpose segmented the image. This takes the longest, expect to see a single period for a small while per image.
1. After the segmented `.seg` file is saved

For the first 5 images and for every 24th image, after running the Cellpose segmenter, the script will tell you about that image. Here's an example:
```
	About image ...[image name]...
		73 cells identified.
		Average Luminance:    Background: 94    Cells: 98
		Area:    Background: 1191228
			 Cells (total): 231524    Cells (average): 3172
	Total time elapsed: 50 seconds
```

All of this information (except time elapsed) is saved to the output file CSV, so you don't have to keep track of the terminal output if you don't want to. 

# Command Line Arguments

To read documentation about the command line arguments, run `python libbylab_cellpose.py -h` with no other command line arguments.

Here's an example of the output: 

```
(cellpose) abigailbishop@vision cellpose % python libbylab_cellpose.py -h

Loading Packages ... Done!
usage: libbylab_cellpose.py [-h] [--cytoplasm_color CYTOPLASM_COLOR]
                            [--nucleus_color NUCLEUS_COLOR]
                            [--cell_diameter CELL_DIAMETER]
                            [--um_per_pixel UM_PER_PIXEL]
                            [--update_csv UPDATE_CSV]
                            [--edge_buffer EDGE_BUFFER]
                            [--image_extension IMAGE_EXTENSION]
                            [--save_boundary_images]
                            dir_to_analyze model_type

positional arguments:
  dir_to_analyze        Path to folder containing images to analyze.
  model_type            Model to use. Give Cellpose model or path to your own.

options:
  -h, --help            show this help message and exit
  --cytoplasm_color CYTOPLASM_COLOR
                        Color of cytoplasm. 0=Gray, 1=Red, 2=Green, 3=Blue
  --nucleus_color NUCLEUS_COLOR
                        Color of nucleus. 0=Gray, 1=Red, 2=Green, 3=Blue
  --cell_diameter CELL_DIAMETER
                        Expected diameter of cells, in pixels.
  --um_per_pixel UM_PER_PIXEL
                        Micrometers per pixel.
  --update_csv UPDATE_CSV
                        If provided, code will update csv output in
                        --dir_to_analyze.
  --edge_buffer EDGE_BUFFER
                        ROIs within this many pixels of the edge are ignored in
                        the area, width, and height calculation.
  --image_extension IMAGE_EXTENSION
                        Extension of the images you intend to analyze.
  --save_boundary_images
                        If present, will save pngs depicting rectangular
                        boundaries used to calculate cell area, width, and height.
```