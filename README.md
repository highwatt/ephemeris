# Ephemeris
Simple C++ header only library for working with EPM and DE ephemeries.

# Usage

```cpp
#include "eph.h"
int main() {

	///test all DE ephemeris
	std::filesystem::path base = "D:\\Projects\\eph\\ascii";
	eph::ephemeris().import_and_test(base / "de102");   
	eph::ephemeris().import_and_test(base / "de200");   
	eph::ephemeris().import_and_test(base / "de202");   
	eph::ephemeris().import_and_test(base / "de403");   
	eph::ephemeris().import_and_test(base / "de405");   
	eph::ephemeris().import_and_test(base / "de406");   
	eph::ephemeris().import_and_test(base / "de410");   
	eph::ephemeris().import_and_test(base / "de421");   
	eph::ephemeris().import_and_test(base / "de422");   
	eph::ephemeris().import_and_test(base / "de423");   
	eph::ephemeris().import_and_test(base / "de424");   
	eph::ephemeris().import_and_test(base / "de430");    
	eph::ephemeris().import_and_test(base / "de430t");   
	eph::ephemeris().import_and_test(base / "de431");   //need ~3.3.Gb memory allocated
	eph::ephemeris().import_and_test(base / "de432");   
	eph::ephemeris().import_and_test(base / "de432t");  
	eph::ephemeris().import_and_test(base / "de433");   
	eph::ephemeris().import_and_test(base / "de434");   
	eph::ephemeris().import_and_test(base / "de435");   
	eph::ephemeris().import_and_test(base / "de436");   
	eph::ephemeris().import_and_test(base / "de436t");  
	eph::ephemeris().import_and_test(base / "de438");   
	eph::ephemeris().import_and_test(base / "de438t");  
	eph::ephemeris().import_and_test(base / "de440");   
	eph::ephemeris().import_and_test(base / "de440t");  
	eph::ephemeris().import_and_test(base / "de441");   //need ~3.5.Gb memory allocated

	/// simple usage
	eph::ephemeris e;
	e.import_and_test(base / "de440t");
	e.save(base / "de440t" / "440t.bin");
	e.load(base / "de440t" / "440t.bin");
	e.test(base / "de440t" / "testpo.440t");

	double x[6] = { 0 };
	e.state(2520576.5, eph::targets::earth, x);

	for (auto val : x)
		std::cout << val << std::endl;

	system("pause");

return 0;
}
```
