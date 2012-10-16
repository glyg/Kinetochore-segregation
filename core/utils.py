import sys

def circ_box(xy, rad):
    return (xy[0] - rad, xy[1] - rad , xy[0] + rad , xy[1] + rad)

def scale(x, size, pix_size = 0.0645):
    """
    Scale the position x on a line of size "size" from microns to pixels
    origin is put in the center of the line
    """
    return int(x / pix_size + size / 2)

def progress(percent):
    """
    Print a progress bar
    percent = -1 to end and remove the progress bar
    """

    # If process is finished
    if percent == -1:
        sys.stdout.write("\r" + " " * 80 + "\r\r")
        sys.stdout.flush()
        return

    # Size of the progress bar
    size = 50

    # Compute current progress
    progress = (percent+1) * size / 100

    # Build progress bar
    bar = "["
    for i in range(progress-1):
        bar += "="
    bar += ">"
    for i in range(size - progress):
        bar += " "
    bar += "]"

    # Write progress bar
    sys.stdout.write("\r%d%% " %(percent+1) + bar)
    sys.stdout.flush()