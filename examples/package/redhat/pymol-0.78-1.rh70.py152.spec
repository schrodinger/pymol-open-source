Summary: PyMOL Molecular Graphics System
Name: pymol
Version: 0.78
Release: 1.rh70.py152
Copyright: Python
Group: Development/Tools
URL: http://www.pymol.org
Source: /usr/src/redhat/SOURCES/pymol-0_78-src.tgz
BuildRoot: /var/tmp/pymol-root

%description
PyMOL, a molecular graphics program with a Python API.

%prep

%setup

%build
make
make pmw
make compileall

%install
rm -rf ${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}/usr
mkdir -p ${RPM_BUILD_ROOT}/usr/bin
mkdir -p ${RPM_BUILD_ROOT}/usr/lib
mkdir -p ${RPM_BUILD_ROOT}/usr/lib/python1.5
mkdir -p ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol
cp -dprv modules ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol/
cp -dprv test ${RPM_BUILD_ROOT}/usr/lib/python1.5/site-packages/pymol/
install -m 755 pymol.com ${RPM_BUILD_ROOT}/usr/bin/pymol

mkdir -p ${RPM_BUILD_ROOT}/usr/doc
mkdir -p ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 LICENSE ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 DEVELOPERS ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 CHANGES ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
install -m 644 sample.pymolrc ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/
cp -dprv examples ${RPM_BUILD_ROOT}/usr/doc/%{name}-%{version}-%{release}/

%clean
rm -rf ${RPM_BUILD_ROOT}

%files
/usr/lib/python1.5/site-packages/pymol/modules/
/usr/lib/python1.5/site-packages/pymol/test/
/usr/doc/%{name}-%{version}-%{release}/LICENSE
/usr/doc/%{name}-%{version}-%{release}/DEVELOPERS
/usr/doc/%{name}-%{version}-%{release}/CHANGES
/usr/doc/%{name}-%{version}-%{release}/sample.pymolrc
/usr/doc/%{name}-%{version}-%{release}/examples/
/usr/bin/pymol


%changelog
* Tue Feb 05 2002 Warren DeLano <warren@delanoscientific.com>
- the clueless leading the blind

