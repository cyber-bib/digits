#include <array>

#include <limits>
#include <type_traits>

#include <exception>
#include <cstdlib>
#include <string>
#include <chrono>
#include <iomanip>
#include <iostream>

template<class T>
struct half_adder {
  // static_assert(std::is_integral_v<T>, "type specified must be an integral type");
  // static_assert(std::is_unsigned_v<T>, "type specified must be an unsigned type");
  // static_assert(  !std::is_const_v<T>, "type specified must not be const type");

  static std::array<T, 2> evaluate(
    T &lhs,
    T &rhs,
    T const &_max = std::numeric_limits<T>::max()
  ) {
    static T const xval = std::numeric_limits<T>::max();
    std::array<T, 2> dst;
    T avg, tmp;

    // lhs = lhs > _max ? (lhs % (_max + 1)): lhs;
    // rhs = rhs > _max ? (rhs % (_max + 1)): rhs;

    avg = lhs/2 + rhs/2;
    avg += (rhs%2 != 0 && lhs%2 != 0 && rhs*lhs != 0) ? 1 : 0;

    tmp = (_max + 1)/2;
    tmp += ((_max + 1)%2 == 0) ? 0 : 1;

    auto q = avg < tmp;

    dst[0] = (q || _max == xval) ? (lhs + rhs) : ((lhs + rhs) % (_max + 1));
    dst[1] = q ? 0 : 1;


    return dst;
  }

  typedef T value_type;

  //inputs
  value_type m_lhs;
  value_type m_rhs;
  value_type m_max;

  //outputs
  value_type m_cout;
  value_type m_sum;

  std::array<value_type, 2> operator()(
    value_type const &lhs,
    value_type const &rhs,
    value_type const &_max = std::numeric_limits<value_type>::max()
  ) {
    this->m_lhs = lhs;
    this->m_rhs = rhs;
    this->m_max = _max;

    auto dst = this->evaluate(
      this->m_lhs,
      this->m_rhs,
      this->m_max  );


    this->m_sum  = dst[0];
    this->m_cout = dst[1];

    return dst;
  }
};

template<class T>
struct full_adder {
  // static_assert(std::is_integral_v<T>, "type specified must be an integral type");
  // static_assert(std::is_unsigned_v<T>, "type specified must be an unsigned type");
  // static_assert(  !std::is_const_v<T>, "type specified must not be const type");

  static std::array<T, 2> evaluate(
    T &lhs,
    T &rhs,
    T &cin = 0,
    T const &_max = std::numeric_limits<T>::max()
  ) {
    std::array<T, 2> dst;

    auto tmp = half_adder<T>::evaluate(lhs, rhs, _max);

    cin = cin > _max ? (cin % (_max + 1)): cin;

    dst = half_adder<T>::evaluate(tmp[0], cin, _max);

    dst[1] += tmp[1];

    return dst;
  }

  typedef T value_type;

  //inputs
  value_type m_lhs;
  value_type m_rhs;
  value_type m_cin;
  value_type m_max;

  //outputs
  value_type m_sum;
  value_type m_cout;

  std::array<value_type, 2> operator()(
    value_type const &lhs,
    value_type const &rhs,
    value_type const &cin = 0,
    value_type const &_max = std::numeric_limits<value_type>::max()
  ) {
    this->m_lhs = lhs;
    this->m_rhs = rhs;
    this->m_cin = cin;
    this->m_max = _max;

    auto dst = this->evaluate(
      this->m_lhs,
      this->m_rhs,
      this->m_cin,
      this->m_max  );


    this->m_sum  = dst[0];
    this->m_cout = dst[1];

    return dst;
  }
};

template<class T>
struct arithmetic_base {
  // deduced virtuals
  // virtual T&   operator=                           ()       {}
  // virtual T    operator+ (arithmetic_base const &rhs) const {}
  // virtual T    operator- (arithmetic_base const &rhs) const {}
  // virtual T    operator~ (arithmetic_base const &rhs) const {}
  // virtual bool operator! ()                           const {}
  // virtual bool operator&&(arithmetic_base const &rhs) const {}
  // virtual bool operator||(arithmetic_base const &rhs) const {}

  // pure virtuals
  #define BUILD_OP_BASE(x, rtype, arg, def) virtual rtype operator##x(arg) def;

  #define BUILD_OP(x) BUILD_OP_BASE(x, T&, T const &rhs, =0)

  BUILD_OP(+= );
  BUILD_OP(-= );
  BUILD_OP(*= );
  BUILD_OP(/= );
  BUILD_OP(%= );
  BUILD_OP(&= );
  BUILD_OP(|= );
  BUILD_OP(^= );
  BUILD_OP(<<=);
  BUILD_OP(>>=);

  #undef BUILD_OP

  #define BUILD_OP(x) BUILD_OP_BASE(x, bool, T const &rhs, const = 0)

  BUILD_OP(< );

  #undef BUILD_OP

  #define BUILD_OP(x) BUILD_OP_BASE(x, T&, void, = 0)

  BUILD_OP(++);
  BUILD_OP(--);

  #undef BUILD_OP
  #undef BUILD_OP_BASE

  // deduced virtuals
  #define BUILD_OP_BASE(x, rtype, arg, qual, code)   \
    virtual rtype operator##x  (arg) qual {          \
      rtype dst = static_cast<rtype const &>(*this); \
      code;                                          \
      return dst;                                    \
    }

  #define BUILD_OP(x) \
    BUILD_OP_BASE(x, T, T const &rhs, const, dst x= rhs)

  BUILD_OP(+ );
  BUILD_OP(- );
  BUILD_OP(* );
  BUILD_OP(/ );
  BUILD_OP(% );
  BUILD_OP(& );
  BUILD_OP(| );
  BUILD_OP(^ );
  BUILD_OP(<<);
  BUILD_OP(>>);

  #undef BUILD_OP

  #define BUILD_OP(x) \
    BUILD_OP_BASE(x, T, int, , x(static_cast<T&>(*this)))

  BUILD_OP(++);
  BUILD_OP(--);

  #undef BUILD_OP
  #undef BUILD_OP_BASE

  #define BUILD_OP(x, code)                          \
    virtual bool operator##x  (T const &rhs) const { \
      bool dst;                                      \
      auto &lhs = static_cast<T const &>(*this);     \
      code;                                          \
      return dst;                                    \
    }

  BUILD_OP(== , dst = !(lhs <  rhs) && !(lhs > rhs));
  BUILD_OP(!= , dst = !(lhs == rhs));
  BUILD_OP(>  , dst =  (rhs <  lhs));
  BUILD_OP(<= , dst = !(lhs >  rhs));
  BUILD_OP(>= , dst = !(lhs <  rhs));

  #undef BUILD_OP

};

template<class T>
struct tiny : arithmetic_base<tiny<T>> {
  T m_val : 5;

  constexpr tiny() : m_val(0) {}
  constexpr tiny(int val) : m_val(val) {}
  tiny& operator=(int val) {
    this->m_val = val;
    return *this;
  }
  friend std::ostream& operator<<(std::ostream& os, tiny<T> const t) {
    os << (short)t.m_val;
    return os;
  }

  #define BUILD_OP_BASE(op, out, in, qual, def) \
    out operator##op(in) qual { def; return *this; };

  #define BUILD_OP(op) \
    BUILD_OP_BASE(op, tiny&, , , op##this->m_val)

  //prefix
  BUILD_OP(++)
  BUILD_OP(--)

  #undef BUILD_OP

  #define BUILD_OP(op) \
    BUILD_OP_BASE(op, tiny&, tiny const &rhs, override, this->m_val op rhs.m_val)

  BUILD_OP(+= )
  BUILD_OP(-= )
  BUILD_OP(*= )
  BUILD_OP(/= )
  BUILD_OP(%= )
  BUILD_OP(&= )
  BUILD_OP(|= )
  BUILD_OP(^= )
  BUILD_OP(<<=)
  BUILD_OP(>>=)

  bool operator<  (tiny const &rhs) const override {
    return this->m_val < rhs.m_val;
  }
};

typedef tiny<char>          stiny;
typedef tiny<unsigned char> utiny;

namespace std {
  template<>
  struct numeric_limits<utiny> {
    static constexpr utiny max() noexcept {
      constexpr utiny dst(-1);
      return dst;
    }
  };
};

bool test_adder() {
  for(int m = 2; m < 32; m++) {
    for(int i = 0; i < 32; i++) {
      for(int j = 0; j < 32; j++) {
        utiny a = i, b = j, s = 0;

        a %= m;
        b %= m;

        auto avg = a/2 + b/2;
        avg += (a%2 != 0 && b%2 != 0 && a*b != 0) ? 1 : 0;

        s = (a+b)%m;

        s = (avg < 16) ? s : s + 32-m;

        int r = s.m_val;

        if(r != ((i+j)%m)) return false;
        // std::cout << std::setw(3) << (i+j)%m;
        // std::cout << std::setw(3) << s;
      }
    }
  }

  return true;
}
int main(int argc, char *argv[]) {
  try {

    // int m = std::stoll(argv[1]);

    std::cout << "test " << (test_adder() ? "passed" : "failed") << "!\n";
    // for(int i = 0; i < 32; i++) {
    //   for(int j = 0; j < 32; j++) {
    //     utiny a = i, b = j, s = 0;

    //     a %= m;
    //     b %= m;

    //     auto avg = a/2 + b/2;
    //     avg += (a%2 != 0 && b%2 != 0 && a*b != 0) ? 1 : 0;

    //     s = (a+b)%m;

    //     s = (avg < 16) ? s : s + 32-m;

    //     std::cout << std::setw(3) << s;
    //   }
    //   std::cout << std::endl;
    // }

    // std::cout << std::endl;

    // for(int i = 0; i < 32; i++) {
    //   for(int j = 0; j < 32; j++) {
    //     // int m = 24;

    //     // auto s = (i + j)%m;
    //     // s = s > (m/2) ? s - m : s;
    //     std::cout << std::setw(3) << (i+j)%m;
    //   }
    //   std::cout << std::endl;
    // }
  }
  catch(...) {
    std::cerr << "Oops!!! something went wrong\n";
    return -1;
  }
  return 0;
}