Please use version inside folder "03-18-2016_release".

Instructions:

1. Place all files inside the "jars" folder inside the following directory: 

Program Files > ImageJ > plugins


2. Simply drag moco_.jar onto ImageJ.

Alternatively (if that doesn't work for some reason),

Open ImageJ, click

Plugins > Install Plugin...

and choose 10092015_release>moco_.jar.


3. Restart ImageJ, and and the moco plugin should show up under Plugins.


---------------------------------
Directions: You must have open the stack and the template image that you want to align the stack with (i.e. the first image in the stack, or the average image of all frames in the stack).



w: The maximimum distance (in pixels) to be translated in the x and y directions. (i.e. w = 2 means that the maximum translation can be 2 pixels up/down
and 2 pixels left/right).

Downsample value: The amount of times to downsample by 2 for faster processing. (i.e. downsampling 512x512 image by 1 makes it 256x256). Leads to problems if downsample value is too large.



Generate log file: Generates a results table that shows the amount by which each image was translated. Can be saved and used for other images/channels.

Plot: Plots the sqrt(x^2 + y^2), i.e. the amount of total translation for each frame
(to plot RMS, you must also choose to generate the log file).


If an image in the stack to be registered is translated by amounts that seem too large, then the program will warn you and will allow you to normalize the intensity in the image (which can improve results).


We recommend that you use 8-bit stacks for registration (there are sometimes issues with 16-bit stacks that lead to translation terms that are too large as noted above). If you generate a log file and save it, you can use it later to register the original 16-bit stack really fast. 


Email guevara.james@gmail.com to ask any questions or report any bugs or feature requests.

