# С++ Эфемериды
[![en](https://img.shields.io/badge/lang-en-blue.svg)](https://github.com/highwatt/ephemeris/blob/main/README.md)

Простая библиотека, написанная на C++ 20 для чтения эфемерид EPM, созданных Институтом прикладной астрономии Российской академии наук (ИПА РАН) и DE эфемерид, созданных JPL NASA. Текущая версия библиотеки была протестирована с эфемеридами **EPM2021**, а также с DE эфемеридами, начиная с **DE102** и заканчивая **DE441**. Описание используемых в библиотеке алгоритмов можно найти в публикации [Практическое использование эфемерид EPM и DE. Труды МАИ, 2022, № 125](https://doi.org/10.34759/trd-2022-125-18).

# Использование эфемерид
Скачайте файлы с данными для DE или EPM эфемерид в формате ASCII:
* [JPL Development Ephemeris](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/ "Jet Propulsion Laboratory (JPL)").
* [Эфемериды ИПА РАН](https://ftp.iaaras.ru/pub/epm/EPM2021/DE/ "Institute of Applied Astronomy of the Russian Academy of Sciences (IAARAS)").

Для работы библиотеки файлы с данными не должны иметь расширений ****.txt*** or ****.bin***. Эти расширения зарезервированы для работы библиотеки. Если скаченные файлы имеют такие раширения их необходимо переименовать в любые другие, но расширения фалов должны совпадать. Например, эфемериды **EPM2021** состоят из двух файлов: ***epm2021_de_ascii.txt*** и ***header.21***. Просто переименуйте файл ***epm2021_de_ascii.txt*** в ***epm2021_de_ascii.21***. Ниже приведен код для импорта данных EPM и DE эфемерид и печати эфемерид Земли на определенную дату.

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

Импорт ASCII файлов может быть очень медленным. Для решения этой проблемы предусмотрена возможность сохранения данных эфемерид в двоичный файл. Загрузка двоичного файла происходит гораздо быстрее и он также имеет меньший размер.
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
Вы также можете использовать для даты и времени синтаксический "сахар" в виде литералов из пространства имен **std::chrono** стандартной библиотеки.
```cpp
	using namespace std::chrono_literals;
	using namespace std::chrono;	
	
	/// Get ephemeris data using std::chrono format of dates & time
	auto [x, y, z, vx, vy, vz] = e.state(eph::date::jd(May/30/2022, 5h+10min), eph::targets::moon);
```
