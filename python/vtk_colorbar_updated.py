import vtk
'''
Helper functions for the creation of colorbar actors in VTK
'''

__all__ = [
    'ColorbarParam',
    'Colorbar',
]

class ColorbarParam:
    def __init__(self, title='No title', title_col=[1,1,1], title_font_size=22, 
        label_col=[1,1,1], pos=[0.9, 0.5], width=80, height=400, 
        nlabels=4, font_size=18, title_offset=10):
        self.title=title
        self.title_col=title_col
        self.label_col=label_col
        self.pos=pos
        self.width=width
        self.height=height
        self.nlabels=nlabels
        self.font_size=font_size
        self.title_offset=title_offset
        self.title_font_size=title_font_size

class Colorbar:
    def __init__(self, ctf, param=None, is_float=True):
        # size and relative position
        if param is None:
            param = ColorbarParam()
        if is_float:
            format='%0.2f'
        else:
            format='%0.0f'  
        # Create a color bar
        self.scalar_bar = vtk.vtkScalarBarActor(position=param.pos, 
                        maximum_width_in_pixels=param.width, 
                        maximum_height_in_pixels=param.height,
                        title=param.title,
                        vertical_title_separation=param.title_offset,
                        number_of_labels=param.nlabels,
                        text_pad=8,
                        draw_tick_labels=True,
                        label_format=format, 
                        unconstrained_font_size=True)
        self.scalar_bar.SetLookupTable(ctf)
        # title properties
        title_prop = vtk.vtkTextProperty(color=param.title_col, shadow=False, 
                        font_size=param.title_font_size, bold=True)
        self.scalar_bar.SetTitleTextProperty(title_prop)
        self.scalar_bar.SetLabelFormat(format)
        label_prop = vtk.vtkTextProperty(color=param.label_col, shadow=False, 
                        font_size=param.font_size, bold=False)
        self.scalar_bar.SetLabelTextProperty(label_prop)

    def get(self):
        return self.scalar_bar
    
    def set_title(self, title, size):
        self.scalar_bar.title = title
        self.scalar_bar.title_text_property.font_size = size

    def set_label(self, nlabels, size):
        self.scalar_bar.number_of_labels = nlabels
        self.scalar_bar.label_text_property.font_size = size

    def set_position(self, pos):
        self.scalar_bar.position = (pos[0], pos[1])

    def set_size(self, width, height):
        self.scalar_bar.maximum_width_in_pixels = width
        self.scalar_bar.maximum_height_in_pixels = height

