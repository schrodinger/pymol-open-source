from enum import IntEnum

from pymol.Qt import QtGui, QtWidgets
from pymol.Qt import QtCore
Qt = QtCore.Qt
QSI = QtGui.QStandardItem  # For brevity


class SceneDictIndex(IntEnum):
    QPIXMAP = 0
    MESSAGE = 1
    ACTIONS = 2


class SceneTableColumn(IntEnum):
    NAME = 0
    IMAGE = 1
    MESSAGE = 2
    ACTIONS = 3


class ScenePanel(QtWidgets.QWidget):
    '''
    Scene Panel dialog for displaying and altering scenes in PyMOL.
    '''

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self)
        self.cmd = parent.cmd
        self.scene_list = self._update_scene_list()
        # EXAMPLE: scene_dict = {name:[QPixmap,Message,Actions]}
        self.scene_dict = {}

        self._build_table_elements(parent)
        self._populate_data()

        self.resize(365, 700)

    def _build_table_elements(self, parent):
        '''
        Create the various elements in the table and add them to the
        proper layout.
        '''
        self.setWindowTitle('Scene Panel')
        layout = QtWidgets.QVBoxLayout(self)
        self.setLayout(layout)

        # Called mid because there will be others
        mid_layout = QtWidgets.QGridLayout()
        layout.addLayout(mid_layout)

        # Mid Elements
        self.sceneTableWidget = QtWidgets.QTableWidget(self)
        mid_layout.addWidget(self.sceneTableWidget, 0, 0)
        self.sceneTableWidget.viewport().installEventFilter(self)

    def eventFilter(self, source, event):
        '''
        Event filter for capturing and processing events.
        '''
        focus_event = 24  # This will be replaced with the event itself.
        if event.type() == focus_event:
            self._populate_data()

        return super().eventFilter(source, event)

    def _populate_data(self):
        '''
        Main loop for populating the table and class data structures.
        '''
        # Clear rows
        self.sceneTableWidget.setRowCount(0)
        self.scene_list = self._update_scene_list()
        self._update_scene_dict()

        self._format_condensed_widget()

        for scene_name in self.scene_list:
            # Add a row for the scene
            row_number = self.sceneTableWidget.rowCount()
            self.sceneTableWidget.insertRow(row_number)
            self.sceneTableWidget.setRowHeight(row_number, 100)
            # Get data from scene_dict
            scene_image_label = self._get_scene_image_label(
                self.scene_dict[scene_name][SceneDictIndex.QPIXMAP])
            scene_message = self.scene_dict[scene_name][SceneDictIndex.MESSAGE]
            scene_actions = self.scene_dict[scene_name][SceneDictIndex.ACTIONS]

            # Set items within the new row
            self.sceneTableWidget.setItem(
                row_number, SceneTableColumn.NAME,
                QtWidgets.QTableWidgetItem(scene_name))
            self.sceneTableWidget.setCellWidget(
                row_number, SceneTableColumn.IMAGE, scene_image_label)
            self.sceneTableWidget.setItem(
                row_number, SceneTableColumn.MESSAGE,
                QtWidgets.QTableWidgetItem(scene_message))
            self.sceneTableWidget.setItem(
                row_number, SceneTableColumn.ACTIONS,
                QtWidgets.QTableWidgetItem(scene_actions))

        self._format_table()

    def _update_scene_dict(self):
        '''
        Updates the class' scene_dict to reflect the current state of
        the class' scene_list. Currently, place holder value are being used
        for the message and action but these will be replaced by functions
        responsible for fetching these values for a given scene.

        Note: This is also currently being used to make sure that scenes
        are being properly removed from the dictionary if they no longer
        appear in the scene_list. There may be a better way of doing this or
        this may be handled in separate step in the future.
        '''
        for scene_position in range(len(self.scene_list)):
            scene_name = self.scene_list[scene_position]
            if scene_name not in self.scene_dict:
                scene_pix_map = self._get_scene_png(scene_name)
                scene_message = 'This is a base message'
                scene_actions = 'Rock, zoom, something'
                self.scene_dict[scene_name] = [
                    scene_pix_map, scene_message,
                    scene_actions, scene_position]

        remove_scenes = [
            scene_name for scene_name in self.scene_dict
            if scene_name not in self.scene_list]

        for scene_name in remove_scenes:
            self.scene_dict.pop(scene_name)

    def _get_scene_png(self, scene_name):
        '''
        Returns a QPixmap object that can be assigned to a label
        in order to display the scene thumbnail.
        This will be obsolete and removed in the future.

        Note:
        Currently, this function goes to the scene and then calls cmd.png to
        generate the PNG buffer.
        This functionality should be transferred to a
        different function that is called when creating a
        new scene from anywhere and the PNG buffer should be
        stored with the scene information. There is a
        subtask set up for this purpose.
        '''
        scene_pix_map = QtGui.QPixmap()
        self.cmd.scene(scene_name, animate=0)
        scene_pix_map.loadFromData(self.cmd.png(None, 200, 100), "PNG")
        self.cmd.refresh()
        return scene_pix_map

    def _get_scene_image_label(self, scene_pix_map):
        '''
        Takes in the QPixmap from get_scene_png and uses it to create a label
        displaying the thumbnail.
        In the future, this function will generate its
        own QPixmap object from the PNG buffer stored with the scene.
        '''
        scene_image_label = QtWidgets.QLabel(self)
        scene_image_label.setPixmap(scene_pix_map)
        return scene_image_label

    def _format_condensed_widget(self):
        '''
        Formats the widget to look and size properly. There will be a
        matching function for the expanded view of the widget.
        '''
        self.sceneTableWidget.setMinimumWidth(275)
        self.setMinimumWidth(375)
        self.resize(300, self.height())
        self.sceneTableWidget.setColumnCount(2)
        self.sceneTableWidget.setHorizontalHeaderLabels(
            ['Name', 'Scene Preview'])
        self.sceneTableWidget.resize(260, 500)

    def _format_table(self):
        '''
        Formats the table to look properly. Currently, vertical headers
        are being hidden but they will be used for reordering.
        '''
        hh = self.sceneTableWidget.horizontalHeader()
        hh.setStretchLastSection(True)
        self.sceneTableWidget.setFocus()
        self.sceneTableWidget.hide()
        self.sceneTableWidget.resizeColumnsToContents()
        self.sceneTableWidget.verticalHeader().setVisible(False)
        self.sceneTableWidget.show()

    def _update_scene_list(self):
        '''
        Returns an updated list of the current scenes. Currently, this
        is very basic but some handling will be done here in the future
        or this will likely be removed.
        '''
        list_scenes = self.cmd.get_scene_list()
        return list_scenes
