'''
Lighting Settings Plugin

(c) 2013-2017 Schrodinger Inc.
'''

#import Tkinter
from pymol import cmd, plugins
from pymol.Qt import QtGui, QtCore, QtWidgets

Qt = QtCore.Qt

class SettingSlider(QtWidgets.QSlider):

    def __init__(self, parent, setting, min_val, max_val, res, line_edit):
        super(SettingSlider, self).__init__(Qt.Horizontal, parent)

        self.setting = setting
        self.min_val = float(min_val)
        self.value_range = float(max_val) - self.min_val
        self.line_edit = line_edit
        self.setObjectName(self.setting)

        self.setMinimum(0)
        self.setMaximum(self.value_range / float(res))

        val = float(cmd.get(setting))
        self.setDoubleValue(val)
        self.updateLineEdit()

        self.valueChanged.connect(self.updateLineEdit)
        self.valueChanged.connect(self.updateSetting)
        self.line_edit.editingFinished.connect(self.lineEditUpdated)

    def lineEditUpdated(self):
        try:
            val = float(self.line_edit.text())
        except:
            self.updateLineEdit()
            return

        self.setDoubleValue(val)
        self.updateLineEdit()

    def update(self, value):
        self.setDoubleValue(value)
        self.updateLineEdit()

    def updateLineEdit(self):
        self.line_edit.setText("%g" % self.getDoubleValue())

    def updateSetting(self):
        cmd.set(self.setting, self.getDoubleValue())

    def getDoubleValue(self):
        return self.min_val + self.value_range * (
            (self.value()) / float(self.maximum()))

    def setDoubleValue(self, val):
        self.setValue(self.maximum() * (val - self.min_val) / self.value_range)


def update_setting(name, value):
    cmd.set(name, value)
    slider = dialog.findChild(SettingSlider, name)
    if slider:
        slider.update(value)

def preset_default():
    update_setting('ambient', 0.14)
    update_setting('direct', 0.45)
    update_setting('spec_direct', 0)
    update_setting('spec_direct_power', 55.)
    update_setting('light_count', 2)
    update_setting('shininess', 55.)
    update_setting('reflect', 0.45)
    update_setting('spec_count', -1)
    update_setting('spec_power', -1.)
    update_setting('spec_reflect', -1.)
    update_setting('specular', 1)
    update_setting('specular_intensity', 0.5)
    update_setting('ambient_occlusion_mode', 0)
    update_setting('ambient_occlusion_scale', 25.0)
    update_setting('ambient_occlusion_smooth', 10)
    update_setting('power', 1.)
    update_setting('reflect_power', 1.)

def preset_metal():
    update_setting('ambient', 0.2)
    update_setting('direct', 0.) # diffuse
    update_setting('spec_direct', 0)
    update_setting('shininess', 51.) # same as spec_power
    update_setting('reflect', 0.5) # diffuse
    update_setting('spec_count', -1)
    update_setting('spec_reflect', -1.)
    update_setting('specular', 1)
    update_setting('specular_intensity', 0.5)  # same as specular

def preset_plastic():
    update_setting('ambient', 0.)
    update_setting('direct', 0.2) # diffuse
    update_setting('spec_direct', 0)
    update_setting('shininess', 32.) # same as spec_power
    update_setting('reflect', 0.55) # diffuse
    update_setting('spec_count', -1)
    update_setting('spec_reflect', -1.)
    update_setting('specular', 1)
    update_setting('specular_intensity', 0.5)  # same as specular

def preset_rubber():
    update_setting('ambient', 0.05)
    update_setting('direct', 0.2) # diffuse
    update_setting('spec_direct', 0)
    update_setting('shininess', 10.) # same as spec_power
    update_setting('reflect', 0.5) # diffuse
    update_setting('spec_count', -1)
    update_setting('spec_reflect', -1.)
    update_setting('specular', 1)
    update_setting('specular_intensity', 0.5)  # same as specular

def preset_xray():
    update_setting('ambient', 1.0)
    update_setting('direct', -0.6)
    update_setting('spec_direct', 0)
    update_setting('spec_direct_power', 55.)
    update_setting('light_count', 2)
    update_setting('shininess', 0.)
    update_setting('reflect', -0.6)
    update_setting('spec_count', -1)
    update_setting('spec_power', -1.)
    update_setting('spec_reflect', -1.)
    update_setting('specular', 0.0)
    update_setting('specular_intensity', 0.0)
    update_setting('ambient_occlusion_mode', 0)
    update_setting('ambient_occlusion_scale', 25.0)
    update_setting('ambient_occlusion_smooth', 10)
    update_setting('power', 1.)
    update_setting('reflect_power', 1.)

dialog = None

def lightingsettings():

    global dialog

    if not dialog:
        dialog = create_dialog()

    dialog.show()
    dialog.raise_()

def create_dialog():

    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle('Lighting Settings')

    setting = plugins.get_pmgapp().skin.setting

    sliders = [
        "Diffuse Reflection",
        ('ambient', 0, 1, None),
        ('reflect', -1, 1, None),

        "Direct Light from Front",
        ('direct (+reflect)', -1, 1, None), # diffuse, coupled with "reflect"
        ('spec_direct', 0, 1, None),
        ('spec_direct_power', 0, 100, 1),

        "Free placeable directed Lights",
        ('light_count', 1, 8, 1),
        ('edit_light', 1, 7, 1),

        "Specular Reflection",
        ('spec_count', -1, 8, 1),
        # ('spec_power', -1, 200, 1), # deprecated since v1.5
        ('shininess', 0, 100, None), # same as spec_power
        ('spec_reflect', -0.01, 1, None),
        ('specular', 0, 1, None),
        ('specular_intensity (=specular)', 0, 1, None), # same as specular

        "Ambient Occlusion (Surface only)",
        ('ambient_occlusion_mode', 0, 2, 1),
        ('ambient_occlusion_scale', 1.0, 50., None),
        ('ambient_occlusion_smooth', 1, 20, 1),

        "Ray trace only",
        ('power', 1, 10, None),
        ('reflect_power', 1, 10, None),
    ]

    layout = QtWidgets.QVBoxLayout(dialog)

    button_layout = QtWidgets.QHBoxLayout()
    layout.addLayout(button_layout)
    layout.setContentsMargins(5, 0, 5, 0)
    button_layout.addWidget(QtWidgets.QLabel("<font color=red>Presets:</font>"))

    presets = [
        ( "Default", preset_default ),
        ( "Metal", preset_metal ),
        ( "Plastic", preset_plastic ),
        ( "Rubber", preset_rubber ),
        ( "X-Ray", preset_xray ),
    ]

    for name, fun in presets:
        btn = QtWidgets.QPushButton(name, dialog)
        btn.pressed.connect(fun)
        btn.setAutoDefault(False)
        button_layout.addWidget(btn)

    form_layout = QtWidgets.QFormLayout()
    form_layout.setContentsMargins(0, 0, 0, 0)
    form_layout.setVerticalSpacing(0)
    form_layout.setLabelAlignment(Qt.AlignLeft)
    layout.addLayout(form_layout)

    for i, item in enumerate(sliders, 1):
        if isinstance(item, str):
            label = QtWidgets.QLabel("<font color=blue>"+item+"</font>")
            form_layout.addRow(label)
            continue

        name, min, max, res = item
        if res is None:
            res = 0.01 if (max - min < 100) else 0.1

        line_edit = QtWidgets.QLineEdit(dialog)
        slider = SettingSlider(dialog,
            name.split()[0], min, max, res, line_edit)

        h_layout = QtWidgets.QHBoxLayout()
        h_layout.addWidget(slider, 3)
        h_layout.addWidget(line_edit, 1)

        form_layout.addRow(name, h_layout)

    return dialog
