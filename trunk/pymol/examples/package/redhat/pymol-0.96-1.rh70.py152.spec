Summary: PyMOL Molecular Graphics System
Name: pymol
Version: 0.96
Release: 1.rh70.py152
Copyright: Python
Group: Development/Tools
URL: http://www.pymol.org
Source: /usr/src/redhat/SOURCES/pymol-0_96-src.tgz
BuildRoot: /var/tmp/pymol-root
Requires: python >= 1.5.2
Requires: tcl >= 8.0.0
Requires: tk >= 8.0.0
Requires: libpng >= 1.0.0
Requires: tkinter >= 1.5.2

%description
PyMOL is a molecular graphics system with an embedded Python
interpreter designed for real-time visualization and rapid generation
of high-quality molecular graphics images and animations. It is fully
extensible and available free to everyone via the "Python"
license. Although a newcomer to the field, PyMOL can already be used
to generate stunning images and animations with ease. It can also
perform many other valuable tasks (such as editing PDB files) to
assist you in your research.

%prep

%setup
cp setup/Rules.linux-rpm-rh70-py152 Rules.make
cp setup/pymol.com.linux-rpm-rh70-py152 pymol.com

%build
cp contrib/modules/ExtensionClass.c-py152 contrib/modules/ExtensionClass.c
cp contrib/modules/ExtensionClass.h-py152 contrib/modules/ExtensionClass.h
make
make pmw
make compileall

%install
rm -rf ${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol
cp -dprv modules ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol/
cp -dprv data ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol/
cp -dprv scripts ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol/
echo pymol/modules > ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol.pth
echo pymol/modules/pymol >> ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol.pth
mkdir -p ${RPM_BUILD_ROOT}/usr/bin
install -m 755 pymol.com ${RPM_BUILD_ROOT}/usr/bin/pymol

mkdir -p ${RPM_BUILD_ROOT}/usr/share/doc/%{name}-%{version}/
install -m 644 LICENSE ${RPM_BUILD_ROOT}/usr/share/doc/%{name}-%{version}/
install -m 644 DEVELOPERS ${RPM_BUILD_ROOT}/usr/share/doc/%{name}-%{version}/
install -m 644 CHANGES ${RPM_BUILD_ROOT}/usr/share/doc/%{name}-%{version}/
install -m 644 setup/sample.pymolrc ${RPM_BUILD_ROOT}/usr/share/doc/%{name}-%{version}/
cp -dprv examples ${RPM_BUILD_ROOT}/usr/share/doc/%{name}-%{version}/

%clean
rm -rf ${RPM_BUILD_ROOT}

%files
/usr/lib/python1.5/site-packages/pymol/
/usr/lib/python1.5/site-packages/pymol.pth
/usr/share/doc/%{name}-%{version}/
/usr/bin/pymol


%changelog

