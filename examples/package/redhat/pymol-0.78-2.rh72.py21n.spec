Summary: PyMOL Molecular Graphics System
Name: pymol
Version: 0.78
Release: 2.rh72.py21n
Copyright: Python
Group: Development/Tools
URL: http://www.pymol.org
Source: /usr/src/redhat/SOURCES/pymol-0_78-src.tgz
BuildRoot: /var/tmp/pymol-root
Requires: python2 >= 2.1
Requires: Numeric >= 20
Requires: tcl >= 8.0.0
Requires: tk >= 8.0.0
Requires: libpng >= 1.0.0

%description
PyMOL is a molecular graphics system with an embedded Python interpreter 
designed for real-time visualization and rapid generation of high-quality 
molecular graphics images and animations. It is fully extensible and 
available free to everyone via the "Python" license. Although a newcomer 
to the field, PyMOL can already be used to generate stunning images and 
animations with unprecedented ease. It can also perform many other 
valuable tasks (such as editing PDB files) to assist you in your research.

%prep

%setup
cp Rules.linux-rpm-rh72-py21n Rules.make
cp pymol.com.linux-rpm-rh72-py21 pymol.com

%build
make
make pmw
make compileall

%install
rm -rf ${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}/usr/lib/python2.1/site-packages/pymol
cp -dprv modules ${RPM_BUILD_ROOT}/usr/lib/python2.1/site-packages/pymol/
cp -dprv test ${RPM_BUILD_ROOT}/usr/lib/python2.1/site-packages/pymol/
cp -dprv examples ${RPM_BUILD_ROOT}/usr/lib/python2.1/site-packages/pymol/
mkdir -p ${RPM_BUILD_ROOT}/usr/bin
install -m 755 pymol.com ${RPM_BUILD_ROOT}/usr/bin/pymol

mkdir -p ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 LICENSE ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 DEVELOPERS ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 CHANGES ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 sample.pymolrc ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/

%clean
rm -rf ${RPM_BUILD_ROOT}

%files
/usr/lib/python2.1/site-packages/pymol/
/usr/doc/%{name}-%{version}-%{release}/LICENSE
/usr/doc/%{name}-%{version}-%{release}/DEVELOPERS
/usr/doc/%{name}-%{version}-%{release}/CHANGES
/usr/doc/%{name}-%{version}-%{release}/sample.pymolrc
/usr/bin/pymol


%changelog
* Wed Feb 06 2002 Warren DeLano <warren@delanoscientific.com>
- modified to use pristine source tarball
* Wed Feb 06 2002 Stephen Adler <adler@bnl.gov>
- cleaned up script and made it python2 complient
* Tue Feb 05 2002 Warren DeLano <warren@delanoscientific.com>
- the clueless leading the blind

