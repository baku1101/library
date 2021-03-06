// auto mod int
// https://youtu.be/L8grWxBlIZ4?t=9858
// https://youtu.be/ERZuLAxZffQ?t=4807 : optimize
// https://youtu.be/8uowVvQ_-Mo?t=1329 : division
// const int MOD = 1000000007;
struct mint {
	LL x; // typedef long long LL;
	mint(LL x=0):x((x%MOD+MOD)%MOD){}
	mint operator-() const { return mint(-x);}
	mint& operator+=(const mint a) {
		if ((x += a.x) >= MOD) x -= MOD;
		return *this;
	}
	mint& operator-=(const mint a) {
		if ((x += MOD-a.x) >= MOD) x -= MOD;
		return *this;
	}
	mint& operator*=(const mint a) {
		(x *= a.x) %= MOD;
		return *this;
	}
	mint operator+(const mint a) const {
		mint res(*this);
		return res+=a;
	}
	mint operator-(const mint a) const {
		mint res(*this);
		return res-=a;
	}
	mint operator*(const mint a) const {
		mint res(*this);
		return res*=a;
	}
	mint pow(LL t) const {
		if (!t) return 1;
		mint a = pow(t>>1);
		a *= a;
		if (t&1) a *= *this;
		return a;
	}

	// for prime MOD
	mint inv() const {
		return pow(MOD-2);
	}
	mint& operator/=(const mint a) {
		return (*this) *= a.inv();
	}
	mint operator/(const mint a) const {
		mint res(*this);
		return res/=a;
	}
};
typedef vector<mint> VM;
typedef vector<VM> VVM;
