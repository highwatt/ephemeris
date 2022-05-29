# Ephemeris
Simple C++ 20 header only library for reading EPM and DE ephemeries. Current version was tested for all DE ephemeris from **DE102** to **DE441** and for **EPM2021** ephemeris.

# Usage
Download ASCII data files for DE or EPM ephemeris:
* [Development Ephemeris by JPL](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/ "Jet Propulsion Laboratory (JPL)").
* [Ephemeris of Planets and Moon by IAARAS](https://ftp.iaaras.ru/pub/epm/EPM2021/DE/ "Institute of Applied Astronomy of the Russian Academy of Sciences (IAARAS)").

ASCII data files must have't ****.txt*** or ****.bin*** extensions. If there are, rename them. File extensions must match. For example **EPM2021** ephemeris have two files: ***epm2021_de_ascii.txt*** and ***header.21***. Just rename ***epm2021_de_ascii.txt*** to ***epm2021_de_ascii.21***. 

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
	/// And test
	e.test(base / "de440t" / "testpo.440t");

	double x[6] = { 0 };
	/// Get ephemeris data for Earth
	e.state(2520576.5, eph::targets::earth, x);

	for (auto val : x)
		std::cout << val << std::endl;

	system("pause");

return 0;
}
```
