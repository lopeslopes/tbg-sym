import numpy as np
from pathlib import Path
import holoviews as hv
from holoviews.operation.datashader import rasterize
from PIL import Image, ImageChops
import logging
logging.getLogger('bokeh').setLevel(logging.FATAL)


hv.extension("bokeh")


def trim(image_path, output_path=None, border_color=(255, 255, 255)):
    img = Image.open(image_path).convert("RGB")
    bg = Image.new("RGB", img.size, border_color)
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox:
        cropped = img.crop(bbox)
        if output_path:
            cropped.save(output_path)
        return cropped
    else:
        print("No borders found (image likely fully solid or transparent).")
        return img


# load data
def load_lattice(path):
    with open(path) as f:
        return np.array([
            [np.float64(x) for x in line.split(";")]
            for line in f if line.strip()
        ])


# load properties
def load_properties(path):
    p = 0
    q = 0
    i = 0
    steps = 0
    with open(path / "properties.dat") as f:
        data = f.readlines()
        for line in data:
            if line != "\n":
                aux = line.split("=")
                if aux[0] == "p":
                    p = int(float(aux[1]))
                elif aux[0] == "q":
                    q = int(float(aux[1]))
                elif aux[0] == "i":
                    i = int(float(aux[1]))
                elif aux[0] == "steps":
                    steps = int(float(aux[1]))
    return p, q, i, steps


def moire_plot(moire_period, folder_name, latA1, latB1, latA2, latB2, latAA, latBA, latAB, latBB):
    r = moire_period
    theta = np.linspace(0, 2 * np.pi, 200)
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    circle_data = np.stack([x, y], axis=1)
    circle_data2 = np.stack([2.0*x, 2.0*y], axis=1)

    # holoviews plot
    # folder_name = str(dataset_dirs[k])
    first_underline = folder_name.find("_")
    second_underline = folder_name.find("_", first_underline+1)
    third_underline = folder_name.find("_", second_underline+1)
    angle_name = folder_name[first_underline+1:second_underline]

    all_pts = np.concatenate((latA1, latB1, latA2, latB2))
    points_total = hv.Points(all_pts)

    min_x = -850
    max_x =  850
    min_y = -850
    max_y =  850
    f_width = max_x-min_x
    f_height = max_y-min_y
    moire = rasterize(points_total, width=f_width+100, height=f_height+100)
    moire = moire.opts(
        cmap="blues",
        # colorbar=True,
        yaxis=None,
        width=1500,
        height=1500,
        xlim=(min_x, max_x),
        ylim=(min_y, max_y),
        title=angle_name,
        aspect="equal"
    )

    overlap = hv.Overlay([])
    if (latAB.size > 0):
        pointsAB = hv.Points(latAB, label="AB").opts(size=15, color="purple")
        overlap *= pointsAB
    if (latBA.size > 0):
        pointsBA = hv.Points(latBA, label="BA").opts(size=15, color="orange")
        overlap *= pointsBA
    if (latAA.size > 0):
        pointsAA = hv.Points(latAA, label="AA").opts(size=15, color="green")
        overlap *= pointsAA
    if (latBB.size > 0):
        pointsBB = hv.Points(latBB, label="BB").opts(size=15, color="magenta")
        overlap *= pointsBB

    circ1 = hv.Curve(circle_data).opts(color="gray")
    circ2 = hv.Curve(circle_data2).opts(color="gray")

    overlap = overlap * circ1 * circ2

    full_image = moire * overlap
    full_image.opts(legend_position="bottom_left")
    hv.save(full_image, "results/500k_"+angle_name+".png", fmt="png", backend="bokeh")
    trim("results/500k_"+angle_name+".png", "results/500k_"+angle_name+".png")


# load available data folders
base_data_path = Path("data")
dataset_dirs = sorted([d for d in base_data_path.iterdir() if d.is_dir()])

latA1 = load_lattice(base_data_path / "latticeA1_500k.dat")
latB1 = load_lattice(base_data_path / "latticeB1_500k.dat")

for k in range(len(dataset_dirs)):
    # obtain lattice properties
    print("Angle ", k, " of ", len(dataset_dirs))
    a = 2.46
    p, q, j, steps = load_properties(dataset_dirs[k])
    angle_i = np.acos((3.0*(q**2) - (p**2))/(3.0*(q**2) + (p**2)))
    angle_f = np.acos((3.0*((q-1)**2) - (p**2))/(3.0*((q-1)**2) + (p**2)))
    angle = angle_i + (j/steps)*(angle_f - angle_i)
    print("Angle in radians: ", angle)
    print("Angle in degrees: ", (angle * 180) / np.pi)
    moire_period = a/(2*np.sin(angle/2))
    print("D = ", moire_period)

    # load data from the lattices

    latA2 = load_lattice(dataset_dirs[k] / "latticeA2.dat")
    latB2 = load_lattice(dataset_dirs[k] / "latticeB2.dat")
    latAA = load_lattice(dataset_dirs[k] / "latticeAA.dat")
    latBA = load_lattice(dataset_dirs[k] / "latticeBA.dat")
    latAB = load_lattice(dataset_dirs[k] / "latticeAB.dat")
    latBB = load_lattice(dataset_dirs[k] / "latticeBB.dat")

    folder_name = str(dataset_dirs[k])
    moire_plot(moire_period, folder_name, latA1, latB1, latA2, latB2, latAA, latBA, latAB, latBB)
