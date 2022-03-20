//---------------------------------------------------------------------------
#pragma once
//---------------------------------------------------------------------------
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cassert>
#include <string_view>
#include <format>
//---------------------------------------------------------------------------
namespace eph {

	struct streambuf_view : public std::streambuf {

		streambuf_view(const char* s, std::size_t count) {
			auto p = const_cast<char*>(s);
			this->setg(p, p, p + count);
		}
		streambuf_view(std::string_view str) :
			streambuf_view(str.data(), str.size()) {
		}
	};	

	/// @brief ������ std::istringstream ������� ����� ������� �� ������ ��� �����������
	class istringstream_view : private virtual streambuf_view, public std::istream {
	public:
		istringstream_view(std::string_view str)
			: streambuf_view(str)
			, std::istream(this) {
		}
	};

	/// @brief ��������� ����. ��� MSVC ���������� ���������� ������. 
	FILE* openfile(std::string_view filename, char const* mode) {

		FILE* file = nullptr;
#ifdef _MSC_VER
		fopen_s(&file, filename.data(), mode);
#else
		file = fopen(file_name.data(), mode);
#endif
		return file;
	}

	/// @brief ������ ���������� ����� � ������
	/// @param filename - ��� ����� ��� ������
	/// @param str - ������ ��� ������
	/// @exeption ������� ���������� std::runtime_error � ������ ������������� ������� ����.
	void getcontent(std::string_view filename, std::string& str) {
		
		FILE* file = openfile(filename, "rb");
		if (!file) throw std::runtime_error(std::string(filename) + " - file not found");

		fseek(file, 0, SEEK_END);
		const size_t size = ftell(file);
		fseek(file, 0, SEEK_SET);

		str.resize(size);
		fread(str.data(), 1, size, file);
		fclose(file);
	}

	/// @brief �������������� �������� ��� ������� ����� ���� ��������� ���������: 
	/// 0 - Mercury, 1 - Venus, 2 - Earth, 3 - Mars, 4 - Jupiter, 
	///	5 - Saturn, 6 - Uranus, 7 - Neptune, 8 - Pluto, 9 - Moon, 10 - Sun, 
	///	11 - Solar-System barycenter, 12 - Earth-Moon barycenter, 17 - Moon (geocentric),
	/// 13 - Earth Nutations, 14 - Lunar mantle libration, 15 - Lunar mantle angular velocity,
	/// 16 - Terrestrial Time (TT) - Barycentric Dynamical Time (TDB), 17 - Moon (geocentric)
	enum class targets {

		mercury,				///  0 - Mercury
		venus,					///  1 - Venus
		earth,					///  2 - Earth
		mars,					///  3 - Mars
		jupiter,				///  4 - Jupiter
		saturn,					///  5 - Saturn
		uranus,					///  6 - Uranus
		neptune,				///  7 - Neptune
		pluto,					///  8 - Pluto
		moon,					///  9 - Moon
		sun,					///  10 - Sun
		solar_system_barycenter,///  11 - Solar-System barycenter
		earth_moon_barycenter,	///  12 - Earth-Moon barycenter
		earth_nutation,			///  13 - Earth Nutations in longitude and obliquity
		lunar_libration,		///  14 - Lunar mantle libration
		lunar_angular_velocity,	///  15 - Lunar mantle angular velocity
		delta_time,				///  16 - Terrestrial Time (TT) - Barycentric Dynamical Time (TDB) at geocenter
		moon_geocentric,		///  17 - Moon (geocentric)
		total				    ///  18
	};

	/// @brief ��������� �������� �������� �� ������� order
	/// @param t - ��������������� ����� [-1, 1]
	/// @param order - ������� �������� ��������
	/// @param pos - ���������, ������������ �������� �������� ��� ���������� ���������
	/// @param vel - ���������, ������������ �������� �������� ��� ���������� �������� 
	inline void fill_poly(double t, size_t order, double pos[], double vel[]) {

		pos[0] = 1.0;
		pos[1] = t;
		vel[0] = 0.0;
		vel[1] = 1.0;
		vel[2] = 4.0 * t;

		for (size_t i = 2; i < order; ++i) {

			pos[i] = 2.0 * t * pos[i - 1] - pos[i - 2];
			vel[i + 1] = 2.0 * t * vel[i] + 2.0 * pos[i] - vel[i - 1];
		}
	}

	/// @brief ��� ������ � ����������� ���������� ������� ������� ���� ephemeris.
	/// ��������: eph::ephemeris e; e.import_and_test("path to directory with ascii data files");
	struct ephemeris {

		static constexpr double epsilon = 1e-13;
		/// @brief ��������� ��������� ��� ������� target.
		/// ��� target = 0 - Mercury, 1 - Venus, 2 - Earth, 3 - Mars, 4 - Jupiter, 
		///	5 - Saturn, 6 - Uranus, 7 - Neptune, 8 - Pluto, 9 - Moon, 10 - Sun, 
		///	11 - Solar-System barycenter, 12 - Earth-Moon barycenter, 17 - Moon (geocentric)
		/// ������ ����������� �������� [x, y, z, vx, vy, vz], ���������� � ��, �������� � ��/�.
		/// 
		/// ��� target = 13 - Earth Nutations
		/// ������ ����������� �������� [dpsi, depsilon, 0, 0, 0, 0], ���.
		/// 
		/// ��� target = 14 - Lunar mantle libration
		/// ������ ����������� �������� [phi, theta, psi, 0, 0, 0], ���.
		/// 
		/// ��� target = 15 - Lunar mantle angular velocity
		/// ������ ����������� �������� [omega_x,omega_y,omega_z, 0, 0, 0], ���/����.
		/// 
		/// ��� target = 16 - Terrestrial Time (TT) - Barycentric Dynamical Time (TDB)
		/// ������ ����������� �������� [dt, 0, 0, 0, 0, 0], �.
		/// 
		/// @param jd - ��������� ����, �� ������� ������������� ���������.
		/// @param target - ������������� ������� ��� �������� ������������� ���������.
		/// @param x - ������ �����������.
		/// @exception ������� ���������� std::out_of_range � ������ target >= 18(total).
		void state(double jd, targets target, double x[6]) {

			switch (target) {
				case targets::mercury:	calculate(jd, series::mercury, x); break;
				case targets::venus:	calculate(jd, series::venus, x); break;
				case targets::earth: {

					double moon_geo[6];
					double factor = 1.0 / (1.0 + constant("EMRAT"));

					calculate(jd, series::earth_moon_barycenter, x);
					calculate(jd, series::moon_geocentric, moon_geo);					

					for (size_t i = 0; i < 6; ++i)
						x[i] -= moon_geo[i] * factor;

					break;
				}
				case targets::mars:		calculate(jd, series::mars, x); break;
				case targets::jupiter:	calculate(jd, series::jupiter, x); break;
				case targets::saturn:	calculate(jd, series::saturn, x); break;
				case targets::uranus:	calculate(jd, series::uranus, x); break;
				case targets::neptune:	calculate(jd, series::neptune, x); break;
				case targets::pluto:	calculate(jd, series::pluto, x); break;
				case targets::moon: {

					double earth[6];
					double factor = 1.0 / (1.0 + constant("EMRAT"));

					calculate(jd, series::earth_moon_barycenter, earth);
					calculate(jd, series::moon_geocentric, x);

					for (size_t i = 0; i < 6; ++i) {
						earth[i] -= x[i] * factor;
						x[i] += earth[i];
					}

					break;
				}
				case targets::sun:						calculate(jd, series::sun, x); break;
				case targets::solar_system_barycenter:	memset(&x[0], 0, 6 * sizeof(double)); break;
				case targets::earth_moon_barycenter:	calculate(jd, series::earth_moon_barycenter, x); break;
				case targets::earth_nutation:			calculate(jd, series::earth_nutation, x); break;
				case targets::lunar_libration:			calculate(jd, series::lunar_libration, x); break;
				case targets::moon_geocentric:			calculate(jd, series::moon_geocentric, x); break;
				case targets::lunar_angular_velocity:	calculate(jd, series::lunar_angular_velocity, x); break;
				case targets::delta_time:				calculate(jd, series::delta_time, x); break;
				default: throw std::out_of_range("Wrong target index."); break;
			}
		}

		/// @brief ������ ��������� ASCII ����� � ����������� �� �������� 
		/// directory � ��������� ���� testpo.xxx
		/// @param directory - ���� � ASCII ������.
		void import_and_test(std::filesystem::path const& directory) {

			import(directory);
			// ���� ���� testpo.xxx
			for (const auto& entry : std::filesystem::directory_iterator(directory)) {
				if (entry.path().stem() == L"testpo") {
					test(entry.path());
					break;
				}
			}
		}

		/// @brief ������ ��������� ASCII ����� � ����������� �� �������� directory
		/// @param directory - ���� � ASCII ������.
		/// @exeption ������� ���������� std::runtime_error � ������ ������������� ������� ����.
		void import(std::filesystem::path const& directory) {

			namespace fs = std::filesystem;

			std::vector<fs::path> file_list;
			fs::path header{ L"header" }, testpo{ L"testpo" }, txt{L".txt"}, bin{L".bin"};

			for (const auto& entry : fs::directory_iterator(directory)) {

				const auto& cur_path = entry.path();
				//��������� ������ ������� ������
				//���������� ����� "header.*", "testpo.*", "*.txt", "*.bin"
				if (cur_path.stem() != header && 
					cur_path.stem() != testpo && 
					cur_path.extension().string() != txt &&
					cur_path.extension().string() != bin)
					file_list.push_back(cur_path);
			}

			std::string ext = file_list.front().extension().string(), content;				
			//��������� ���� ���������				
			std::ifstream file(directory / header.concat(ext), std::ios::binary);				
			//��������� ����� ��������� ����� � ����� ���������� ���������� ��������
			if (!file) {
					
				std::wstring_view add_ext[] = { {L"_571"}, {L"_572"}, {L"_576"} };				
				for (auto& e : add_ext) {
					fs::path new_header = header;
					file.open(directory / new_header += e, std::ios::binary);
					if (file) break;
				}					
				if (!file) throw std::runtime_error(header.string() + " - file not found");
			}

			//������ ���������� ����� � ������
			std::getline(file, content, '\0');
			file.close();

			//������ header.XXX
			read_header(content);

			//��������� ������ ������� �������������
			size_t size = columns_ * (date_to_index(jd_end_) + 1);

			//��� ������� ������� ��������� �������� de431 ��� de441 ��������� �������� ������ 2�� ������
			//��� �������� ������ ��� ���������� x64 ������. ��� x32 ������� ��������������. 
			//����� coefficients_.resize ������ std::bad_alloc. 
			if (size * sizeof(double) > INT_MAX && sizeof(void*) != 8)
				std::cout << std::format("Target platform is x86. Ephemeris de{} require more then 2Gb memory allocation, and can't be read.\n", ext);

			coefficients_.resize(size);

			fact_jd_begin_ = jd_end_, fact_jd_end_ = jd_begin_;

			for (const auto& entry : file_list) {

				std::cout << "Reading: " << entry << std::endl;

				//������ ���������� ����� � ������
				getcontent(entry.string(), content);				

				//�������������� ������ 'D' �� 'e'
				fix_exponent(content);

				//��������� ������ � �����
				istringstream_view stream(content);

				size_t n, local_section;	double begin, end;

				//������ ����� ������ � ���������� �����. � ������
				while (stream >> local_section >> n) {

					if (n != columns_) throw std::out_of_range(entry.string() + " reading error - NCOEFF mismatch");

					//������ ���� ������ � ����� ������
					stream >> begin >> end;

					if (begin < fact_jd_begin_) fact_jd_begin_ = begin;
					if (end > fact_jd_end_) fact_jd_end_ = end;

					size_t section = date_to_index(begin);
					size_t j = section * n, e = j + n;

					//������ ������������
					for (; j < e; ++j) {
						char next = stream.peek();
						//������������ ��������, ����� � ���������� ������ ������ ������ 3-� �������������
						if (next == '\n' && (j % n > n-3))
							break;
						stream >> coefficients_[j];
					}
				}
			}

			if (fact_jd_begin_ != jd_begin_ || fact_jd_end_ != jd_end_)
				std::cout << "WARNING : Missmatch range that ephemeris covers:\nfrom header.xxx Start Epoch : " << 
				std::format("JED = {} Final Epoch : JED = {}\nfrom data files Start Epoch: JED={} Final Epoch: JED={}\n", jd_begin_, jd_end_, fact_jd_begin_, fact_jd_end_);
		}

		/// @brief ��������� ���� �� ����� testpo.xxx
		/// @param file_name - ����, ������� ��� ����� testpo.xxx
		/// @exeption ������� ���������� std::runtime_error � ������ ������������� ������� ����.
		void test(std::filesystem::path const& file_name) {

			std::cout << std::format("Ephemeris test: {}\n", file_name.string());
	
			std::string content;
			//������ ���������� ����� � ������
			getcontent(file_name.string(), content);

			//���� ������ �������� ������
			size_t pos = content.find("EOT");

			std::istringstream stream(content.substr(pos + 3, std::string::npos));

			double jd, test_value;
			double target_data[6], center_data[6];

			size_t n, target, center, index, iter = 0;
			long year, month, day;
			char c;

			std::cout << "--MSG--   de# ---date---  ---jed--- t# c# x# --coordinate--  --difference--\n";

			while (stream >> n >> year >> c >> month >> c >> day >> jd >> target >> center >> index >> test_value) {

				std::string date = std::format("{}.{}.{}", year, month, day);

				if (jd < fact_jd_begin_ || jd > fact_jd_end_) {
					std::cout << std::format("WARNING : {} {} {} {} {} {} {} - Coordinates not found\n", n, date, jd, target, center, index, test_value);
					continue;
				}
				//������������ ������� �.�. � ����� ���������� � �������
				target--, index--;

				//��������� ��� target
				state(jd, targets(target), target_data);

				//��������� ��� center
				if (center != 0) state(jd, targets(--center), center_data);
				else center_data[index] = 0.0;

				double delta = target_data[index] - center_data[index];

				//���� target �� ������� ��� ��������, ��� �������� �������� ������, ��� TT-TDB ��������� � �.�.
				if (target != to_underlying(targets::earth_nutation) && 
					target != to_underlying(targets::lunar_libration) && 
					target != to_underlying(targets::lunar_angular_velocity) && 
					target != to_underlying(targets::delta_time))
					delta /= constant("AU");

				delta = fabs(delta - test_value);

				//���� target - �������� ���������� �������������� psi
				if (target == to_underlying(targets::lunar_libration)) {

					//������������ psi
					if (index == 2) delta /= 0.23 * (jd - 2451545.0);
					//������������  dpsi/dt
					if (index == 5) delta *= 0.01 / (1.0 + (jd - 2451545.0) / 365.25e+2);
				}

				if (delta > epsilon) {
					std::cout << std::format("ERROR   : {} {} {} {} {} {} {} delta={} - Invalid ephemeris\n", n, date, jd, target, center, index, test_value, delta);
					return;
				}
				//����� �� ������ ������ 100 ��������
				if (iter++ % 100 == 0)
					std::cout << std::format("SUCCESS : {} {} {} {} {} {} {} delta={}", n, date, jd, target, center, index, test_value, delta) << std::endl;
			}
			std::cout << std::format("Test '{}' successfully completed.\n\n", file_name.string());
		}

		//������ ��������� �� ��������� �����
		void load(std::filesystem::path const& file_name) {

			FILE* file = openfile(file_name.string(), "rb");
			if (!file) throw std::runtime_error(file_name.string() + " - file not found");

			//���� ������ � ����� ��������
			fread(&jd_begin_, sizeof(jd_begin_), 1, file);
			fread(&jd_end_, sizeof(jd_end_), 1, file);
			fread(&fact_jd_begin_, sizeof(fact_jd_begin_), 1, file);
			fread(&fact_jd_end_, sizeof(fact_jd_end_), 1, file);
			//���������� ���� � ������
			fread(&section_days_, sizeof(section_days_), 1, file);
			//���������� �������� � ������� ������������� ��������� ��������
			fread(&columns_, sizeof(columns_), 1, file);
			fread(&num_targets_, sizeof(num_targets_), 1, file);
			
			size_t n = 0;
			//���������� ��������
			fread(&n, sizeof(n), 1, file);
			//�������� ������ ��� ��������� � �� �����
			const_names_.resize(n);
			const_values_.resize(n);
			//������ ���� ��������
			fread(&const_names_[0], sizeof(Fixbuf), n, file);
			//������ ��������
			fread(&const_values_[0], sizeof(double), n, file);

			//���������� ������������� ��������� ��������
			fread(&n, sizeof(n), 1, file);
			coefficients_.resize(n);
			//������ ������������� ��������� ��������
			fread(&coefficients_[0], sizeof(double), n, file);
			//��������� ��� ������������
			fread(&layout_[0], sizeof(layout_), 1, file);
			fclose(file);
		}

		//��������� ��������� � �������� ����
		void save(std::filesystem::path const& file_name) {

			FILE* file = openfile(file_name.string(), "wb");
			if (!file) throw std::runtime_error(file_name.string() + " - file not found");

			//���� ������ � ����� ��������
			fwrite(&jd_begin_, sizeof(jd_begin_), 1, file);
			fwrite(&jd_end_, sizeof(jd_end_), 1, file);
			fwrite(&fact_jd_begin_, sizeof(fact_jd_begin_), 1, file);
			fwrite(&fact_jd_end_, sizeof(fact_jd_end_), 1, file);
			//���������� ���� � ������
			fwrite(&section_days_, sizeof(section_days_), 1, file);
			//���������� �������� � ������� ������������� ��������� ��������
			fwrite(&columns_, sizeof(columns_), 1, file);
			fwrite(&num_targets_, sizeof(num_targets_), 1, file);

			size_t n = const_names_.size();
			//���������� ��������
			fwrite(&n, sizeof(n), 1, file);
			//������ ���� ��������
			fwrite(&const_names_[0], sizeof(Fixbuf), n, file);
			//������ ��������
			fwrite(&const_values_[0], sizeof(double), n, file);

			n = coefficients_.size();
			//���������� ������������� ��������� ��������
			fwrite(&n, sizeof(n), 1, file);
			//������ ������������� ��������� ��������
			fwrite(&coefficients_[0], sizeof(double), n, file);
			//��������� ��� ������������
			fwrite(&layout_[0], sizeof(layout_), 1, file);
			fclose(file);
		}	

		//�������� ����� ��������� �� �����
		double constant(std::string_view name) const {

			for (size_t i = 0, e = const_names_.size(); i < e; ++i)
				if (name == const_names_[i].buf) return const_values_[i];

			throw std::runtime_error("Constant not found.");
		}

	private:

		enum class series : size_t {
			mercury,				//0 - Mercury			
			venus,					//1 - Venus
			earth_moon_barycenter,	//2 - Earth-Moon barycenter
			mars,					//3 - Mars,
			jupiter,				//4 - Jupiter,
			saturn,					//5 - Saturn,
			uranus,					//6 - Uranus,
			neptune,				//7 - Neptune,
			pluto,					//8 - Pluto,
			moon_geocentric,		//9 - Moon (geocentric)
			sun,					//10 - Sun
			earth_nutation,			//11 - 1980 IAU nutation angles
			lunar_libration,		//12 - Lunar mantle libration (Euler) angles
			lunar_angular_velocity,	//13 - Lunar mantle angular velocity
			delta_time,				//14 - TT - TDB(at geocenter)
			total			
		};	

		template< class Enum >
		static constexpr std::underlying_type_t<Enum> to_underlying(Enum e) noexcept {
			return static_cast<std::underlying_type_t<Enum>>(e);
		}

		enum class poly : size_t  {
			coeff_begin,
			order,
			sub_intervals,
			total
		};		

		void calculate(double jd, series target, double x[6]) {

			if (to_underlying(target) >= num_targets_) throw std::out_of_range("Wrong target index.");

			//�������� ������������� ��� ���� target �� ������ ������
			size_t coeff_begin = layout_[to_underlying(poly::coeff_begin)][to_underlying(target)] - layout_[0][0];

			//������� ��������
			size_t poly_order = layout_[to_underlying(poly::order)][to_underlying(target)];

			//���� � ��������� �� ������� ���������� ��������� �������� ��� (��������)
			//������ �� ������������ �� ����� ������� �������������. time_cover � ���� ���������� �������������   
			size_t time_cover = layout_[to_underlying(poly::sub_intervals)][to_underlying(target)];

			//������ ������, ������� ����������� jd
			//�������� �������, ����� ������ ������ �������������� ��������� (jd1, jd2]
			//� ��������� ������, ��� jd==jd2 date_to_index ������ ������ ��������� ������,
			//��� ��� ��������� ������ �������� � ������ ������� �� ��������
			size_t section = date_to_index(jd - 0.5);

			//���� ������ ������
			double t_section = jd_begin_ + section * section_days_;

			//��������� �������� ������������
			size_t date_subinterval = size_t(section_days_ / time_cover);

			//������ - �������� �� ������ ������
			size_t sub_shift = size_t((jd - t_section - 0.5) / date_subinterval);

			//���� ������ ������������
			double t_sub = t_section + sub_shift * date_subinterval;

			//����������� �������
			size_t dim = target == series::earth_nutation ? 2 : (target == series::delta_time ? 1 : 3);

			//������ ������ �������������
			size_t i = section * columns_ + coeff_begin + sub_shift * dim * poly_order;

			//��������������� ����� [-1, 1]
			double normal_time = 2.0 * (jd - t_sub) / date_subinterval - 1.0;

			//�������� ��� ��������� � ��������
			double pos_poly[max_poly_order], vel_poly[max_poly_order];

			assert(poly_order < max_poly_order && "Invalid Chebyshev polynom order");

			//��������� ��������
			fill_poly(normal_time, poly_order, pos_poly, vel_poly);

			//������� ������
			memset(x, 0, 6 * sizeof(double));

			double* y = &coefficients_[i];

			//���� �� ������������� � �������� �������,
			//�.�. ��� ���������� �������� ���������� �������
			//��������� � ����������� ����������
			for (size_t j = 0; j < dim; ++j) {
				for (size_t k = poly_order; k --> 0; ) {
					//���������
					x[j] += y[j * poly_order + k] * pos_poly[k];
					//��������
					x[j + dim] += y[j * poly_order + k] * vel_poly[k];
				}

				x[j + dim] *= 2.0 / date_subinterval;
			}
		}

		/// @brief ���������� ������ ������, ��������������� ���������� ����.
		/// @param jd - ��������� ����.
		/// @return ������ � ������� coefficients_.
		size_t date_to_index(double jd) {

			return size_t( (jd - jd_begin_) / section_days_ );
		}

		/// @brief ������ ������������ ���� �������� header.XXX.
		/// @param content - ���������� ����� header.XXX.
		void read_header(std::string_view const& content) {		

			size_t pos1010 = content.find("GROUP   1010", 0);
			size_t pos1030 = content.find("GROUP   1030", pos1010);
			size_t pos1040 = content.find("GROUP   1040", pos1030);
			size_t pos1041 = content.find("GROUP   1041", pos1040);
			size_t pos1050 = content.find("GROUP   1050", pos1041);
			size_t pos1070 = content.find("GROUP   1070", pos1050);
			
			//���������� �������� � ����� ������
			size_t n = 12;
			std::string_view head_line = content.substr(0, pos1010);
			std::string_view group1010 = content.substr(pos1010 + n, pos1030 - pos1010 - n);
			std::string_view group1030 = content.substr(pos1030 + n, pos1040 - pos1030 - n);
			std::string_view group1040 = content.substr(pos1040 + n, pos1041 - pos1040 - n);
			std::string_view group1041 = content.substr(pos1041 + n, pos1050 - pos1041 - n);
			std::string_view group1050 = content.substr(pos1050 + n, pos1070 - pos1050 - n);
			
			{//������ ���������� ������������� � ������
				char ncoeff[] = "NCOEFF=";
				size_t pos = head_line.find("NCOEFF=", 0);
				istringstream_view stream(head_line.substr(pos + sizeof(ncoeff)));
				stream >> columns_;
			}
			
			{//������ ������ 1030 ���� ������ � ����� ������
				istringstream_view stream(group1030);
				stream >> jd_begin_ >> jd_end_ >> section_days_;
			}

			{//������ ������ 1040 ����� ��������
				istringstream_view stream(group1040);
				stream >> n;
				const_names_.resize(n);
				for (size_t i = 0; i < n; ++i)
					stream >> const_names_[i].buf;
			}			

			{//������ ������ 1041 �������� ��������, �������������� ������ 'D' �� 'e'
				std::string tmp(group1041);
				fix_exponent(tmp);
				istringstream_view stream(tmp);
				stream >> n;
				const_values_.resize(n);
				for (size_t i = 0; i < n; ++i)
					stream >> const_values_[i];
			}
			
			auto skip_spaces = [](std::istream& st) {
				while (1) {
					char next = st.peek();
					if (next != ' ') return;
					st.get();
				}
			};

			{//������ ������ 1050 ��������� ��� ������������
				istringstream_view stream(group1050);
				// num_targets_ == 0 ��. ���������� num_targets_
				// ������ ������ ������ ������ 1050 �� ������� '\n' � ���������� ����������� �������� num_targets_
				do {
					stream >> layout_[0][num_targets_++];
					//� ��������� ����������, �������� � 410, � ����� ������ ���� �������, �� ����������
					skip_spaces(stream);
				} while (stream.peek() != '\n');

				for (size_t i = 1; i < 3; ++i)
					for (size_t j = 0; j < num_targets_; ++j)
						stream >> layout_[i][j];
			}
		}

		/// @brief ������ 'D'(Fortran) �� 'e'(C++) � ����������� ��c�������.
		/// @param content - ������ � �������.
		void fix_exponent(std::string& content) {
			for (size_t i = 0, n = content.size(); i < n; ++i)
				if (content[i] == 'D') content[i] = 'e';
		}

		//����� ��� �������� ���� ��������
		struct Fixbuf { char buf[7] = { 0 }; };

		//������������ ������� ��������
		static constexpr size_t max_poly_order = 18;
		
		//������ ���� ��������
		std::vector<Fixbuf> const_names_;

		//������ ��������
		std::vector<double> const_values_ = { 0.0 };

		//��������� ������ ������������� ��������� ��������
		//�������� [n x columns_], ��� n - ���������� ������
		//������ ������, ��������������� ���������� ����
		//���������� ������� date_to_index
		std::vector<double> coefficients_ = { 0.0 };

		//��������� ��� ������� � ������������� ��������� �������� ��� ������ �������
		//������� 3�15
		size_t layout_[to_underlying(poly::total)][to_underlying(series::total)] = { 0 };

		//���� ������ � ����� ��������
		double jd_begin_ = { 0.0 }, jd_end_ = { 0.0 }, fact_jd_begin_= { 0.0 }, fact_jd_end_= { 0.0 };

		//���������� ���� � ������
		double section_days_ = { 0.0 };		

		//���������� �������� � ������� ������������� ��������� ��������
		//columns = NCOEF
		size_t columns_ = { 0 };		
		
		// ���������� �������� ��� �������� ����� ���� 13 ��� 15 == poly::total
		// num_targets_ - ����������� ���������� ��������, ������������ �� header.XXX
		size_t num_targets_ = { 0 };
	};	
}//end of epf namespace
//---------------------------------------------------------------------------