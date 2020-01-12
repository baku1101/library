auto isOk = [&](LL mid){
	// ここになにか書く
};

// ng:false，ok:true
LL ng = -1, ok = INT32_MAX;
while (abs(ng - ok) > 1) {
	LL mid = (ok + ng) / 2;
	(isOk(mid) ? ok : ng) = mid;
}

cout << ok << endl;
