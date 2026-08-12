// Copyright (C) 2026 ORIC Labs. GPL-3.0-only.
//
// The underlying FFT is R Core GPL code with process-global scratch metadata.
// Serialize calls so the safe Rust wrapper cannot race that state.
use std::sync::Mutex;

static FFT_LOCK: Mutex<()> = Mutex::new(());

extern "C" {
    fn fft_factor(n: i32, maxf: *mut i32, maxp: *mut i32);
    fn fft_work(
        a: *mut f64,
        b: *mut f64,
        nseg: i32,
        n: i32,
        nspn: i32,
        isn: i32,
        work: *mut f64,
        iwork: *mut i32,
    ) -> bool;
}

pub(crate) fn transform(real: &mut [f64], imaginary: &mut [f64], inverse: bool) {
    assert_eq!(real.len(), imaginary.len());
    let n = i32::try_from(real.len()).expect("FFT length exceeds i32");
    let _guard = FFT_LOCK.lock().expect("R FFT lock poisoned");
    let mut maxf = 0;
    let mut maxp = 0;
    unsafe { fft_factor(n, &mut maxf, &mut maxp) };
    assert!(maxf > 0, "R FFT could not factor transform length {n}");
    let mut work = vec![0.0; 4 * maxf as usize];
    let mut iwork = vec![0; maxp as usize];
    // R stores complex vectors interleaved and calls fft_work with an indexing
    // increment of two. Use the same layout, not merely an equivalent pair of
    // arrays, so the operation order is identical down to floating-point bits.
    let mut interleaved = Vec::with_capacity(real.len() * 2);
    for (&re, &im) in real.iter().zip(imaginary.iter()) {
        interleaved.push(re);
        interleaved.push(im);
    }
    let direction = if inverse { 2 } else { -2 };
    let ok = unsafe {
        fft_work(
            interleaved.as_mut_ptr(),
            interleaved.as_mut_ptr().add(1),
            1,
            n,
            1,
            direction,
            work.as_mut_ptr(),
            iwork.as_mut_ptr(),
        )
    };
    assert!(ok, "R FFT rejected transform length {n}");
    for (index, pair) in interleaved.chunks_exact(2).enumerate() {
        real[index] = pair[0];
        imaginary[index] = pair[1];
    }
}

#[cfg(test)]
mod tests {
    use super::transform;

    #[test]
    fn matches_public_r_fft_observation_and_inverse_scaling() {
        let original = [1.0, 2.0, 3.0, 4.0];
        let mut real = original;
        let mut imaginary = [0.0; 4];
        transform(&mut real, &mut imaginary, false);
        assert_eq!(real, [10.0, -2.0, -2.0, -2.0]);
        assert_eq!(imaginary, [0.0, 2.0, 0.0, -2.0]);

        transform(&mut real, &mut imaginary, true);
        for (observed, expected) in real.into_iter().zip(original) {
            assert!((observed - 4.0 * expected).abs() < 1e-12);
        }
        assert!(imaginary.into_iter().all(|value| value.abs() < 1e-12));
    }
}
