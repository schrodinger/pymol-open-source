Summary: PyMOL Molecular Graphics System
Name: pymol
Version: 0.78
Release: 1.rh70.py152
Copyright: Python
Group: Development/Tools
URL: http://www.pymol.org
Source: /usr/src/redhat/SOURCES/pymol-0_78-src.tgz
BuildRoot: /usr/lib/python1.5/site-packages/pymol

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
cp -dprv modules ${RPM_BUILD_ROOT}/modules
cp -dprv test ${RPM_BUILD_ROOT}/test
install -m 755 pymol.com ${RPM_BUILD_ROOT}/
mkdir -p /usr/doc/%{name}-%{version}-%{release}/
install -m 644 LICENSE /usr/doc/%{name}-%{version}-%{release}/
install -m 644 DEVELOPERS /usr/doc/%{name}-%{version}-%{release}/
install -m 644 CHANGES /usr/doc/%{name}-%{version}-%{release}/
install -m 644 sample.pymolrc /usr/doc/%{name}-%{version}-%{release}/
cp -dprv examples /usr/doc/%{name}-%{version}-%{release}/

rm -rf /usr/local/bin/pymol
ln -s ${RPM_BUILD_ROOT}/pymol.com /usr/local/bin/pymol

%clean
rm -rf ${RPM_BUILD_ROOT}

%files
/modules/
%doc LICENSE
%doc DEVELOPERS
%doc CHANGES
%doc sample.pymolrc
%doc examples/
/pymol.com

%changelog
* Tue Feb 05 2002 Warren DeLano <warren@delanoscientific.com>
- the clueless leading the blind

