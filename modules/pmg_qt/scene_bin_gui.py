from enum import IntEnum
from dataclasses import dataclass

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


@dataclass
class SceneWithTablePosition:
    name: str
    position: int


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
        self._set_button_connections()

        self.resize(365, 700)

    def _build_table_elements(self, parent):
        '''
        Create the various elements in the table and add them to the
        proper layout.
        '''
        self.setWindowTitle('Scene Panel')
        layout = QtWidgets.QVBoxLayout(self)
        self.setLayout(layout)

        top_layout = QtWidgets.QGridLayout()
        layout.addLayout(top_layout)
        mid_layout = QtWidgets.QGridLayout()
        layout.addLayout(mid_layout)
        low_layout = QtWidgets.QGridLayout()
        layout.addLayout(low_layout)

        # Top Elements
        self.instructionLabel = QtWidgets.QLabel(self)
        self.instructionLabel.setText(
            'Double click selected thumbnail to \nload into Workspace.')
        top_layout.addWidget(self.instructionLabel, 0, 0)

        self.addSceneButton = QtWidgets.QPushButton(self)
        self.addSceneButton.setText('Add Scene')
        top_layout.addWidget(self.addSceneButton, 0, 1)

        # Mid Elements
        self.sceneTableWidget = QtWidgets.QTableWidget(self)
        mid_layout.addWidget(self.sceneTableWidget, 0, 0)
        self.sceneTableWidget.viewport().installEventFilter(self)

        self.sceneTableWidget.itemChanged.connect(self._item_changed)
        self.sceneTableWidget.selectionModel().selectionChanged.connect(
            self._selection_changed)
        self.sceneTableWidget.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectRows)
        self.sceneTableWidget.verticalHeader().setSectionsMovable(True)

        # Lower Buttom Elements
        self.deleteButton = QtWidgets.QPushButton(self)
        self.deleteButton.setText("Delete Scene")
        low_layout.addWidget(self.deleteButton, 0, 1)
        self.deleteButton.setEnabled(False)

        self.updateButton = QtWidgets.QPushButton(self)
        self.updateButton.setText("Update Scene")
        low_layout.addWidget(self.updateButton, 0, 0)
        self.updateButton.setEnabled(False)

    def _set_button_connections(self):
        self.addSceneButton.clicked.connect(self._add_scene)
        self.deleteButton.clicked.connect(self._delete_scene)
        self.updateButton.clicked.connect(self._update_scene)
        self.sceneTableWidget.doubleClicked.connect(self._show_selected_scene)

    def eventFilter(self, source, event):
        '''
        Event filter for capturing and processing events.
        '''
        focus_event = 24  # This will be replaced with the event itself.
        paint_event = 12
        if event.type() == paint_event:
            self._check_table_state()
        elif event.type() == focus_event:
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

        self.sceneTableWidget.selectionModel().clearSelection()

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
        '''
        scene_pix_map = QtGui.QPixmap()
        png_buf = self.cmd.get_scene_thumbnail(scene_name)
        scene_pix_map.loadFromData(png_buf, "PNG")
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
        self._set_vertical_headers()
        self.sceneTableWidget.show()

    def _set_vertical_headers(self):
        # Get row count
        row_count = self.sceneTableWidget.rowCount()

        # Iterate in range, setting item
        for ind in range(row_count):
            # Create  arrows item
            arrows = QtWidgets.QTableWidgetItem('\u2195')
            arrows_font = arrows.font()
            arrows_font.setPointSize(20)
            arrows.setFont(arrows_font)
            self.sceneTableWidget.setVerticalHeaderItem(ind, arrows)

    def _update_scene_list(self):
        '''
        Returns an updated list of the current scenes. Currently, this
        is very basic but some handling will be done here in the future
        or this will likely be removed.
        '''
        list_scenes = self.cmd.get_scene_list()
        return list_scenes

    def _add_scene(self):
        '''
        Calls cmd.scene with append and repopulates the table.
        Table is scrolled to the bottom to show new scene.
        '''
        self.cmd.scene('new', 'append', quiet=0)
        self._populate_data()
        self.sceneTableWidget.scrollToBottom()

    def _update_scene(self):
        '''
        Calls cmd.scene with update, replacing the selected scene with
        new scene. Image is not currently being updated but it should.
        Button only active with selection.
        '''
        selection = self.sceneTableWidget.selectionModel().selectedIndexes()
        if selection:
            name = selection[0].data()
            self.cmd.scene(name, 'update')
            self.scene_dict[name][SceneDictIndex.QPIXMAP] = self._get_scene_png(
                name)
            self._update_table()

    def _delete_scene(self):
        '''
        Calls cmd.scene with clear, removing the selected scene.
        Button only active with selection.
        '''
        selection = self.sceneTableWidget.selectionModel().selectedIndexes()
        name = selection[0].data()
        self.cmd.scene(name, 'clear')
        try:
            self.scene_dict.pop(name)
            self.scene_list.remove(name)
        except Exception as e:
            print("Item not found")
            print(e)
        self._update_table()

    def _show_selected_scene(self):
        '''
        Calls cmd.scene with recall, recalling scene to workspace.
        '''
        selection = self.sceneTableWidget.selectionModel().selectedIndexes()
        if selection:
            name = selection[0].data()
            self.cmd.scene(name, 'recall')

    def _update_table(self):
        '''
        Currently pointless but scrolling/formatting will happen
        here to make updating the table more seemless after
        intially populating.
        '''
        self._populate_data()

    def _compare_scene_lists(self, old_list, new_list):
        '''
        Compares current and new list of scenes.
        Returns list of elements that differ ([old, new])
        '''
        diff_list = []
        for i in range(min(len(old_list), len(new_list))):
            if new_list[i] != old_list[i]:
                diff_list.append([old_list[i], new_list[i]])
        return diff_list

    def _get_table_scene_list(self):
        '''
        Returns a sorted list of scene names based on the current
        state of the table. This will be a critical part of reordering,
        but also used for detecting renaming.
        '''
        scene_coor_list = []
        for row_num in range(self.sceneTableWidget.rowCount()):
            name = self.sceneTableWidget.model().data(
                self.sceneTableWidget.model().index(row_num, 0))
            position = int(
                self.sceneTableWidget.rowViewportPosition(row_num)/100)

            scene_coor_list.append(SceneWithTablePosition(name, position))

        scene_coor_list.sort(key=lambda scene: scene.position)

        return [scene_with_pos.name for scene_with_pos in scene_coor_list]

    def _check_table_state(self):
        '''
        Checks the current state of the table for reordering.
        '''
        scene_ordered_list = self._get_table_scene_list()
        self.scene_list = self.cmd.get_scene_list()
        diff_list = self._compare_scene_lists(
            self.scene_list, scene_ordered_list)
        if len(diff_list) > 1:
            # This is the reordering case since multiple have changed
            self._reorder_scenes(scene_ordered_list)

    def _selection_changed(self):
        '''
        Enables buttons when selection is made.
        '''
        item_selected = bool(
            self.sceneTableWidget.selectionModel().selectedIndexes())
        self.deleteButton.setEnabled(item_selected)
        self.updateButton.setEnabled(item_selected)

    def _rename_scene(self, item):
        '''
        Takes in an item from _item_changed and checks conditions before
        calling scene with 'rename'.
        '''
        if ' ' in item.text():
            print("Scene names with spaces are not supported")
            self._update_table()
        elif not item.text():
            print("Blank scene names are not allowed")
            self._update_table()
        else:
            diff_list = self._compare_scene_lists(
                self.cmd.get_scene_list(),
                self._get_table_scene_list())
            if diff_list:
                for old_scene, new_scene in diff_list:
                    self.cmd.scene(old_scene, 'rename', new_key=new_scene)

    def _reorder_scenes(self, scene_ordered_list):
        '''
        Calls cmd.scene_order which requires a string of scene names.
        This string is created using the ordered scene list generated
        by _check_table_state.
        '''
        scene_order_string = ' '.join(
            scene for scene in scene_ordered_list if scene in self.scene_list)
        self.cmd.scene_order(scene_order_string)
        self.scene_list = self._update_scene_list()

    def _item_changed(self, item):
        """
        Called every time an item in the table is changed.
        @param item: The item which has changed
        @type  item: QStandardItem
        """
        if item.column() == SceneTableColumn.NAME:
            self._rename_scene(item)
