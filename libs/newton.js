const version = "1.0.1"; // Oleg Mojivic

const Norm = (v) => {
  return Math.sqrt(
    v.reduce((s, e) => {
      return s + e * e;
    }, 0)
  );
};

const Dot = (a, b) => {
  let s = 0;
  for (const i in a) s += a[i] * b[i];
  return s;
};

// QR - decomposition A = QR of matrix A
const QRDec = (A) => {
  const m = A.length;
  const R = [];

  for (const j in A) {
    R[j] = [];
    for (const i in A) R[j][i] = 0;
  }

  const Q = [];
  for (const i in A) {
    Q[i] = [];
    for (const j in A[0]) Q[i][j] = A[i][j];
  }

  // Q is a copy of A
  for (let i = 0; i < m; i++) {
    let e = Q[i];
    let r = Math.sqrt(Dot(e, e));
    if (r == 0) throw new Error("QRDec: singular matrix");
    R[i][i] = r;
    for (let k in e) e[k] /= r;
    // normalization
    for (let j = i + 1; j < m; j++) {
      const q = Q[j];
      const s = Dot(e, q);
      for (let k in q) q[k] -= s * e[k];
      // orthogonalization
      R[j][i] = s;
    }
  }
  return [Q, R];
};

// QR - backsubstitution
// input: matrices Q,R, array b; output: array x such that QRx=b
const QRBack = (Q, R, b) => {
  const m = Q.length;
  const c = new Array(m);
  const x = new Array(m);
  for (const i in Q) {
    // c = QˆT b
    c[i] = 0;
    for (const k in b) c[i] += Q[i][k] * b[k];
  }
  for (let i = m - 1; i >= 0; i--) {
    // back substitution
    let s = 0;
    for (let k = i + 1; k < m; k++) s += R[k][i] * x[k];
    x[i] = (c[i] - s) / R[i][i];
  }
  return x;
};

// calculates inverse of matrix A
const Inverse = (A) => {
  const t = QRDec(A);
  const Q = t[0];
  const R = t[1];
  const m = [];
  for (const i in A) {
    n = [];
    for (const k in A) {
      n[k] = k == i ? 1 : 0;
    }
    m[i] = QRBack(Q, R, n);
  }
  return m;
};

const Zero = (fs, x, opts = {}) => {
  // Newton's root-finding method

  if (opts.acc == undefined) opts.acc = 1e-6;
  if (opts.dx == undefined) opts.dx = 1e-3;
  if (opts.max == undefined) opts.max = 40; // max iterations
  const J = [];
  for (const i in x) {
    J[i] = [];
    for (const j in x) J[i][j] = 0;
  }

  let minusfx = [];
  let v = fs(x);
  if (v == null) throw "unable to compute fs at " + JSON.stringify(x);
  for (const i in x) minusfx[i] = -v[i];
  do {
    if (opts.max-- < 0) return null;
    for (const i in x)
      for (const k in x) {
        // calculate Jacobian
        x[k] += opts.dx;
        v = fs(x);
        if (v == null) throw "unable to compute fs at " + JSON.stringify(x);
        J[k][i] = (v[i] + minusfx[i]) / opts.dx;
        x[k] -= opts.dx;
      }
    let t = QRDec(J);
    let Q = t[0];
    let R = t[1];
    let Dx = QRBack(Q, R, minusfx);
    // Newton's step
    let s = 2;
    let z = [];
    let minusfz = [];
    do {
      // simple backtracking line search
      s = s / 2;
      z = [];
      for (const i in x) {
        z[i] = x[i] + s * Dx[i];
      }

      minusfz = [];
      v = fs(z);
      if (v == null) throw "unable to compute fs at " + JSON.stringify(z);
      for (const i in x) {
        minusfz[i] = -v[i];
      }
    } while (Norm(minusfz) > (1 - s / 2) * Norm(minusfx) && s > 1 / 128);
    minusfx = minusfz;
    x = z;
    // step done
  } while (Norm(minusfx) > opts.acc);

  return x;
};

const Solve = (
  fs,
  res,
  opts = {
    acc: 1e-6,
    dx: 1e-3,
    max: 40,
  }
) => {
  if (typeof res != "object") res = [typeof res == "number" ? +res : 0];
  const _fs = fs;

  fs = function (x) {
    const r = _fs(x);
    if (typeof r == "number") r = [r];
    for (const i in r) r[i] -= res[i];
    return r;
  };

  const start = [];
  if (opts.start) {
    start = opts.start;
  } else {
    for (const i in res) start[i] = 0;
  }

  const n = Zero(fs, start, opts);
  if (n && n.length == 1) n = n[0];
  return n;
};

const Newton = {
  version,
  Norm,
  Dot,
  QRDec,
  QRBack,
  Inverse,
  Zero,
  Solve,
};
export default Newton;
