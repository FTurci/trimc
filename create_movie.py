import ovito
from ovito.io import *
from ovito.vis import *
import sys

filename = sys.argv[1]
node = import_file(filename, multiple_frames = True, columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z"])
node.add_to_scene()

vp = ovito.dataset.viewports.active_vp

settings = RenderSettings(
    filename = filename+".gif",
    size = (640, 480),
    range = RenderSettings.Range.ANIMATION

)

vp.zoom_all()
vp.render(settings)