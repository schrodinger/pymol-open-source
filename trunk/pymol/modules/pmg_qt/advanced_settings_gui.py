from pymol import setting
from pymol.Qt import QtGui, QtWidgets
from pymol.Qt import QtCore, QtCoreModels
Qt = QtCore.Qt
QSI = QtGui.QStandardItem  # For brevity


class PyMOLAdvancedSettings(QtWidgets.QWidget):
    """
    The Advanced Settings dialog for PyMOL.  This iterates through all possible
    options and adds them to a filterable table.
    """
    def __init__(self, parent, cmd):
        QtWidgets.QWidget.__init__(self, parent, Qt.Window)
        self.setMinimumSize(400, 500)
        self.cmd = cmd

        self.model = QtGui.QStandardItemModel(self)
        self.proxy_model = QtCoreModels.QSortFilterProxyModel(self)
        self.proxy_model.setSourceModel(self.model)

        self.setWindowTitle('PyMOL Advanced Settings')
        layout = QtWidgets.QVBoxLayout(self)
        self.setLayout(layout)
        self.filter_le = QtWidgets.QLineEdit(self)
        layout.addWidget(self.filter_le)
        self.filter_le.setPlaceholderText("Filter")
        self.filter_le.textChanged.connect(self.proxy_model.setFilterRegExp)

        self.populateData()

        self.table = QtWidgets.QTableView(self)
        self.table.setModel(self.proxy_model)
        layout.addWidget(self.table)

        self.formatTable()

        self.model.itemChanged.connect(self.itemChanged)

    def populateData(self):
        """
        Fill the model with data from PyMOL
        """
        from textwrap import fill

        name_list = sorted(setting.get_name_list())

        for name in name_list:
            index = setting._get_index(name)
            v_type, v_list = self.cmd.get_setting_tuple(index)

            value_item = QSI()
            name_item = QSI(name)
            name_item.setFlags(Qt.ItemIsEnabled)
            if v_type == 1:  # CheckBox type
                value_item.setCheckable(True)
                value_item.setEditable(False)  # Can't edit text (but toggles)
                if v_list[0]:
                    value_item.setCheckState(Qt.Checked)
            else:  # Text type
                if v_type in (2, 6):  # int, str
                    value_item.setText(str(v_list[0]))
                elif v_type in (3, 4, 5):  # float, 3f, color
                    value_item.setText(self.cmd.get(index))

            self.model.appendRow([name_item, value_item])
            value_item.setData(index)

    def formatTable(self):
        """
        Set up the table to look appropriately
        """
        hh = self.table.horizontalHeader()
        hh.setStretchLastSection(True)
        hh.setVisible(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setFocus()
        self.table.hide()
        self.table.resizeColumnsToContents()
        self.table.show()

    def itemChanged(self, item):
        """
        Called every time an item in the table is changed, only value items
        are changeable.
        @param item: The item which has changed
        @type  item: QStandardItem
        """

        try:
            index = item.data().toInt()[0]
        except AttributeError:  # PySide not PyQT, different QVariant handling
            index = item.data()

        if item.isCheckable():
            checked = item.checkState() == Qt.Checked
            self.cmd.set(index, checked, log=1, quiet=0)
        else:
            self.cmd.set(index, item.text(), log=1, quiet=0)
