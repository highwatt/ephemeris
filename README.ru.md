# �++ ���������
[![en](https://img.shields.io/badge/lang-en-blue.svg)](https://github.com/highwatt/ephemeris/blob/main/README.md)

������� ����������, ���������� �� C++ 20 ��� ������ �������� EPM, ��������� ���������� ���������� ���������� ���������� �������� ���� (��� ���) � DE ��������, ��������� JPL NASA. ������� ������ ���������� ���� �������������� � ����������� **EPM2021**, � ����� � DE �����������, ������� � **DE102** � ���������� **DE441**. �������� ������������ � ���������� ���������� ����� ����� � ���������� [������������ ������������� �������� EPM � DE. ����� ���, 2022, � 125](https://doi.org/10.34759/trd-2022-125-18).

# ������������� ��������
�������� ����� � ������� ��� DE ��� EPM �������� � ������� ASCII:
* [JPL Development Ephemeris](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/ "Jet Propulsion Laboratory (JPL)").
* [��������� ��� ���](https://ftp.iaaras.ru/pub/epm/EPM2021/DE/ "Institute of Applied Astronomy of the Russian Academy of Sciences (IAARAS)").

��� ������ ���������� ����� � ������� �� ������ ����� ���������� ****.txt*** or ****.bin***. ��� ���������� ��������������� ��� ������ ����������. ���� ��������� ����� ����� ����� ��������� �� ���������� ������������� � ����� ������, �� ���������� ����� ������ ���������. ��������, ��������� **EPM2021** ������� �� ���� ������: ***epm2021_de_ascii.txt*** � ***header.21***. ������ ������������ ���� ***epm2021_de_ascii.txt*** � ***epm2021_de_ascii.21***. ���� �������� ��� ��� ������� ������ EPM � DE �������� � ������ �������� ����� �� ������������ ����.

```cpp
#include "eph.h"
int main() {

	/// Set your path to downloaded ASCII data files 
	std::filesystem::path base = "D:\\Projects\\eph\\ascii";

	double x[6] = { 0 };
	/// Test EPM ephemeris
	eph::ephemeris epm;
	/// Import ASCII files from epm2021 folder
	epm.import(base / "epm2021");
	/// Get ephemeris data for Earth
	epm.state(eph::date::jd(2022, 5, 30), eph::targets::earth, x);

	for (auto val : x)
		std::cout << val << std::endl; 
	
	double y[6] = { 0 };
	/// Test DE ephemeris
	eph::ephemeris de;
	/// Import ASCII files from de440t folder
	de.import(base / "de440t");
	/// Get ephemeris data for Earth
	de.state(eph::date::jd(2022, 5, 30), eph::targets::earth, y);

	for (auto val : y)
		std::cout << val << std::endl;

	system("pause");

return 0;
}
```

������ ASCII ������ ����� ���� ����� ���������. ��� ������� ���� �������� ������������� ����������� ���������� ������ �������� � �������� ����. �������� ��������� ����� ���������� ������� ������� � �� ����� ����� ������� ������.
```cpp
#include "eph.h"
int main() {

	/// Set your path to downloaded ASCII data files 
	std::filesystem::path base = "D:\\Projects\\eph\\ascii";	
	/// Create ephemeris object
	eph::ephemeris e;
	/// Load from ASCII files and test
	e.import_and_test(base / "de440t");
	/// Save to binary file
	e.save(base / "de440t" / "440t.bin");
	/// Load from binary file
	e.load(base / "de440t" / "440t.bin");

	double state[6] = { 0 };
	/// Get ephemeris data for Earth-Moon Barycenter
	e.state(eph::date::jd(2022, 5, 30), eph::targets::earth_moon_barycenter, state);
	
	for (auto val : state)
		std::cout << val << std::endl;
	
	/// Get ephemeris data using structured binding
	auto [x, y, z, vx, vy, vz] = e.state(eph::date::jd(2022, 5, 30), eph::targets::moon);

	std::cout 
		<< "x = " << x << std::endl
		<< "y = " << y << std::endl
		<< "z = " << z << std::endl
		<< "vx = " << vx << std::endl
		<< "vy = " << vy << std::endl
		<< "vz = " << vz << std::endl;	

	system("pause");

return 0;
}
```
�� ����� ������ ������������ ��� ���� � ������� �������������� "�����" � ���� ��������� �� ������������ ���� **std::chrono** ����������� ����������.
```cpp
	using namespace std::chrono_literals;
	using namespace std::chrono;	
	
	/// Get ephemeris data using std::chrono format of dates & time
	auto [x, y, z, vx, vy, vz] = e.state(eph::date::jd(May/30/2022, 5h+10min), eph::targets::moon);
```
