# Detector arrays for a full-lattice heterodyne DM gauge: web survey

**Question.** A heterodyne (beat-frequency, phase-shifting) interferometer must read a
96 x 96 actuator lattice at >= 1 element per actuator (>= ~100 x 100 over the pupil,
ideally 128 x 128 to 512 x 512), each element demodulating a ~1-100 kHz beat
CONTINUOUSLY, to picometer-class phase per second.  What array technologies exist
or are credible within 5-10 years?

**Date of survey:** 2026-10-11.  Web survey only; nothing here was measured.  Every
number carries a source URL; "n/s" = the source does not state it; "n/c" = I could
not confirm it from a primary source.  Numbers marked **[derived]** are my
arithmetic from the sourced inputs, with the formula shown, and are not vendor
claims.

---

## 0. The yardstick used throughout

For a heterodyne phase read in one element, the shot-noise-limited phase
uncertainty is of order `sigma_phi ~ 1/sqrt(N)` rad for `N` detected
photoelectrons in the averaging time (fringe visibility 1, LO >> signal; a
factor ~sqrt(2) either way depending on the LO/signal split -- order of magnitude
only).  At lambda = 633 nm, 1 pm of SURFACE in a double-pass interferometer is
2 pm of OPD = `2*pi*2e-12/633e-9 = 2.0e-5 rad`.  So **"1 pm per second per
element" needs N ~ 2.5e9 photoelectrons per second per element** (1 pm) or
~2.5e7 (10 pm) **[derived]**.  The per-element photon THROUGHPUT -- full well
multiplied by how often the well is emptied, or the count-rate ceiling -- is the
first-order fitness test for every technology below.  Second-order: the readout
quantization / electronics noise must sit below that shot noise, and the whole
array must be read continuously.

Two implications worth stating up front:

* A conventional integrating camera is NOT disqualified by full well alone --
  full well x frame rate is the throughput.  A 20 ke- well emptied 20,000 times
  a second is 4e8 e-/s, pm-class **[derived]**.  What disqualifies most cameras
  is that they cannot stream 4+ samples per beat period over >= 1e4 pixels
  continuously, and that read noise is paid on every sample.
* Any in-pixel lock-in (analog integrator or digital up/down counter) removes
  the per-sample read-noise penalty and the data-rate penalty by integrating
  the demodulated signal for the whole averaging time inside the pixel.  That
  is why class 1 below ranks first.

---

## 1. Lock-in / demodulating CMOS cameras

### 1a. Heliotis heliCam C4 (current product, in-pixel analog I/Q lock-in)

| Item | Value | Source |
|---|---|---|
| Pixels | 542 x 512 with in-pixel lock-in amplifiers (277,504 dual-phase demodulators) | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Pitch / fill factor | 24 um x 24 um; 21 % (S40), 55 % with microlens (S40U/S41U) | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Full well | 500 ke- (C4.0-S40, -S40U); 125 ke- (C4.0-S41U) | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf); [heliotis.com/en/c4](https://www.heliotis.com/en/c4/) |
| Demodulation frequency | 1 kHz - 250 kHz (datasheet); **305 Hz - 134 kHz on the vendor web page** -- discrepancy, see flag | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf); [heliotis.com/en/c4](https://www.heliotis.com/en/c4/) |
| Filter time constant | 0.2 ms - 100 ms | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Internal sample rate | 1 MSa/s per channel (277 GSa/s total) | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Max OUTPUT sample rate | 5.6 kSa/s per channel (3.12 GSa/s total); memory depth 1000 time-resolved output frames | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Output word | 2 x 10 bit (X, Y); outputs magnitude & phase, I & Q, grey | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Phase resolution | 0.1 deg (vendor figure, per sample) | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Background (common-mode) suppression | 21 bit equivalent; min detectable contrast 1e-4 time-resolved, 1e-5 with averaging | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| QE | 61.7 % @400 nm ... 78.7 % @500 nm ... 69.4 % @700 nm, 45.5 % @800 nm | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Read noise | n/s | -- |
| Interface | Gigabit Ethernet; GenICam producer DLL; Python/Matlab/LabVIEW/C++/.NET examples | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| External reference | **"reference input, internal reference output, configurable output, acquisition trigger with programmable phase delay"** -- the demodulation reference CAN be fed from the beat source | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Power / size | 24 V DC, 30 W; 46 x 65 x 180 mm; 500 g | [C4 datasheet](https://www.findlight.net/uploads/products/datasheet/helicam-c4_1675885002.pdf) |
| Price | n/s (quote only) | [findlight product page](https://www.findlight.net/imaging-and-vision/cameras/imaging-spectroscopy/helicam-c4) |
| Announced | 2023 (inVISION news item 2023-07-03) | [inVISION](https://invision-news.com/optische-messtechnik/lock-in-camera-with-demodulation-frequency-of-250khz/) |

**C4M variant (larger format):** 1024 x 1102 pixels, demodulation 305 Hz - 50 kHz,
18-bit-equivalent common-mode suppression, FWC 140 ke- (C4M.2-S4M1U)
([heliotis.com/en/c4](https://www.heliotis.com/en/c4/)).  This is the only
~1 Mpixel lock-in imager I found anywhere.

**Fitness [derived]:**
* Pixel count: 542 x 512 (C4) and 1024 x 1102 (C4M) cover 128 x 128 to 512 x 512
  with margin.  OK.
* Beat frequency: 1-250 kHz (C4) covers the 1-100 kHz band; C4M caps at 50 kHz.  OK.
* Continuous: the in-pixel integrator runs continuously with tc up to 100 ms; the
  I/Q output is the integrated result, so the averaging happens in the pixel.
  At full frame the 5.6 kSa/s figure is a burst-to-memory rate (1000 frames):
  542 x 512 x 20 bit x 5.6 kSa/s = 31 Gbit/s, far above GigE.  Continuous
  streaming over GigE at full frame is ~1e9 / (542 x 512 x 20) ~ 180 I/Q frames/s;
  a 128 x 128 ROI streams at ~3 kHz **[derived]**.  For pm-per-second work you want
  tc ~ 10-100 ms anyway, i.e. 10-100 I/Q frames/s -- GigE is not the limit.
* Photon budget: the datasheet does not say whether 500 ke- is per internal
  1-MSa/s sample or per output sample.  Bounding **[derived]**: per output sample
  at 5.6 kSa/s -> 2.8e9 e-/s -> sigma_phi ~ 1.9e-5 rad/sqrt(s) ~ **1 pm surface
  per sqrt(s)**; per internal sample -> 5e11 e-/s -> 1.4e-6 rad -> 0.07 pm.
  Either bound is pm-class.
* The likely real limit is the **10-bit X/Y output**: 0.1 deg = 1.7e-3 rad per
  output sample.  Averaging 5600 samples/s with dithered quantization gives
  ~2.3e-5 rad/sqrt(s) ~ 1.2 pm surface **[derived]** -- right at the line, and
  only if the quantization noise averages.  Whether the in-pixel integrator can
  be run at tc = 100 ms WITHOUT saturating the 10-bit output on a bright beat
  (i.e. whether gain is adjustable) is n/s and must be asked of the vendor.
* Verdict: **the only commercial device that meets every requirement on paper
  today.**  Key limit: output quantization / unknown read noise; single vendor.

### 1b. Heliotis heliCam C3 (previous generation, still sold)

| Item | Value | Source |
|---|---|---|
| Pixels / pitch | 300 x 300; 40 um | [Liu et al. 2016, PMC4957977](https://pmc.ncbi.nlm.nih.gov/articles/PMC4957977) |
| Demodulation frequency | up to 250 kHz | [Heliotis Swissphotonics presentation](https://www.swissphotonics.net/libraries.files/Heliotis_presentation.pdf) |
| I/Q frame rate | up to 3800 frames/s (paper); "up to 5000 I/Q frames per second for 300x300 pixel" (presentation) | [PMC4957977](https://pmc.ncbi.nlm.nih.gov/articles/PMC4957977); [presentation](https://www.swissphotonics.net/libraries.files/Heliotis_presentation.pdf) |
| ADC | 10-bit | [PMC4957977](https://pmc.ncbi.nlm.nih.gov/articles/PMC4957977) |
| Min integration used | 286 us (user's minimum) | [PMC4957977](https://pmc.ncbi.nlm.nih.gov/articles/PMC4957977) |
| Full well, read noise | n/s | -- |
| Interface | MVTec HALCON interface exists (vendor doc 404 at survey time) | [MVTec](https://www.mvtec.com/products/interfaces/documentation/view/1303-camera-13-mvtecdoc-helicamc3) |

Verdict: 300 x 300 covers the lattice at ~3 px/actuator; superseded by C4.

**Credibility 5-10 y (Heliotis):** four generations shipped (C2 -> C3 -> C4 -> C4M;
the Swissphotonics deck already describes the in-pixel architecture and a
"1 Mpixel -> ~1000 frames/s" direction, which C4M now realizes).  Vendor keeps a
public list of papers using the camera
([heliotis.com/en/helicam_papers](https://www.heliotis.com/en/helicam_papers/)) and
a research page for C4 ([heliotis.com/en/research_c4](https://www.heliotis.com/en/research_c4/));
paper counts per year n/c.  Small Swiss company (Lucerne) -- single-source risk
is the main credibility caveat, not the technology.

### 1c. Other "lock-in" imagers found

* **Optical lock-in camera for GW interferometers (Pockels-cell gating + sCMOS),
  Cao et al. 2019** (arXiv:1907.05224): a Pockels cell + QWP + PBS form an
  optical switch that intensity-modulates the beam at the beat frequency so each
  pixel of an ordinary Andor Zyla 4.2 sCMOS (2048 x 2048, 16-bit, 89 dB dynamic
  range) acts as a lock-in; "2 Mpx resolution at 10 Hz with a sensitivity of
  -62 dBc when averaged over 2s"; demodulated a 36 MHz heterodyne beat in the
  design and 15.4 MHz in the test, four frames per I/Q set
  ([arXiv 1907.05224](https://arxiv.org/abs/1907.05224)).  Lesson: ANY camera can
  be made a lock-in by optical gating, but the photon throughput is still full
  well x frame rate of that camera, and each I/Q needs four frames.
* **Research lock-in pixels** (Kawahito, Toyohashi: lateral-electric-field charge
  modulators) exist in the literature as range-imaging / time-resolved pixels
  ([Kawahito, IISW 2013](https://imagesensors.org/Past%20Workshops/2013%20Workshop/2013%20Papers/10-6_112_Kawahito.pdf));
  no product found.
* No second commercial vendor of a kHz-class in-pixel lock-in imager was found.
  Flag: this is a one-vendor technology class.

---

## 2. Indirect-ToF demodulation sensors (2-tap / N-phase pixels)

These demodulate at MHz, not kHz.  To use one, the heterodyne beat would have to
be moved to the 4-100 MHz band (two AOMs), where the GW phase-camera work shows
a ToF chip does read optical beat phase per pixel.

| Sensor | Pixels / pitch | Modulation range | Taps / phases | Frame rate | Full well | External clock / lock | Status | Sources |
|---|---|---|---|---|---|---|---|---|
| Sony IMX556PLR | 640 x 480 active; 10 um | 4 - 100 MHz (sony-semicon ToF page) | 2 taps (A/B); 4 phases/frame; A-B, A+B, A, B, A&B readout modes | 60 fps full res (4 phase/frame) | "Saturation signal min. 1023 LSB" -- electrons n/s | sensor is the master: "4-wire SPI Master for Illumination signal"; MIPI CSI-2 4-lane 960 Mbps/lane | shipping (industrial) | [flyer PDF](https://www.sony-semicon.com/files/62/flyer_industry/IMX556PLR_PLR1_Flyer.pdf); [Sony iToF page](https://www.sony-semicon.com/en/products/is/industry/tof/itof.html) |
| Melexis MLX75027 | 640 x 480; 10 um CAPD | 4 - 100 MHz in 1 MHz steps, PLL from a fixed 8 MHz external crystal clock | 2-tap; 4 raw phases per distance frame, up to 8 | 120 fps full res (4 phases, Tint 300 us, 4-lane) | single-tap FWC 106 (min) / 160 (typ) ke-; pixel dark noise 83 e-; demod contrast 85 % @40 MHz, 78 % @100 MHz | internal PLL locked to 8 MHz input; sensor drives LEDP/LEDN; MIPI CSI-2 | automotive, shipping | [MLX75027 datasheet (prelim v0.9)](https://www1.futureelectronics.com/doc/Melexis/MLX75027RTC-ABA-200-TR.pdf); [Melexis product page](https://www.melexis.com/en/product/MLX75027/Automotive-VGA-Time-Of-Flight-Sensor) |
| Analog Devices ADSD3100 | 1024 x 1024; 3.5 um; BSI CW-iToF | **n/c** (datasheet fetch timed out twice) | n/c | n/c | n/s | SPI / I2C / MIPI CSI-2 up to 1.5 Gbps/lane | shipping; ADTF3175 module $197 @1k (2022) | [gophotonics summary](https://www.gophotonics.com/products/cmos-image-sensors/analog-devices-inc/21-82-adsd3100); [ADI datasheet PDF](https://www.analog.com/media/en/technical-documentation/data-sheets/adsd3100.pdf); [ISW blog](http://image-sensors-world.blogspot.com/2022/07/high-resolution-tof-module-from-analog.html) |
| Infineon REAL3 IRS2877C / IRS2381C | 640 x 480 (IRS2877C); 224 x 172, 14 um (IRS2381C) | family "up to 100 MHz" | n/s | n/s | n/s | n/s | consumer, shipping | [IRS2877A brief](https://www.infineon.com/dgdl/Infineon-IRS2877AS_REAL3_image_sensor-ProductBrief-v01_00-EN.pdf?fileId=8ac78c8c80f4d32901814971dc03260c); [IRS2381C brief](https://www.mouser.com/pdfDocs/Infineon-REAL3_image_sensor_IRS2381C-PB-v01_00-EN.pdf) |
| pmd flexx2 (dev kit) | 224 x 172 | n/s (two-frequency modes) | n/s | up to 60 fps | n/s | n/s | dev kit | [pmdtec](https://pmdtec.com/en/company/news/announcing-the-next-generation-3d-tof-depth-sensing-development-kit-flexx2/) |
| Teledyne e2v Hydra3D(+) | 832 x 600; 10 um; 3-tap | n/s ("transfer time from 20 ns") | 3 taps | n/s | n/s | n/s | shipping | [Teledyne](https://www.teledynevisionsolutions.com/products/hydra3d_plus/?model=EV3S0M5B-CLVH000-T&vertical=tvs-e2v&segment=tvs) |
| ESPROS epc660 | 320 x 240 | 0.75 - 24 MHz | n/s | n/s | **8,000 ke-** ("ambient and signal") | "external LED/LD modulation input MODCLK ... up to 96 MHz"; 4 MHz system clock internal or external -- whether MODCLK drives the PIXEL demodulation is n/c | shipping | [epc660 datasheet](https://www.neumueller.com/datenblatt/epc/EPC_660.pdf) |

**Metrology precedent:** the Nikhef/VU "high frame-rate phase camera" for
gravitational-wave wavefront sensing uses a TI OPT8241 ToF sensor (320 x 240),
images beats at 9 MHz and 45 MHz in real time at up to 150 fps, reaches
"-62 dBc/second/pixel" single-frame and "-67 dBc/pixel" with 50-frame averaging,
and is "predominantly limited by the electronic noise", shot noise "a factor of
8 below the total noise" ([arXiv 2105.07886](https://arxiv.org/html/2105.07886v1)).
That is the one direct demonstration of a ToF chip as a heterodyne phase
imager, and it was electronics-noise-limited, not photon-limited.

**Fitness [derived]:**
* Pixel count: 0.3-1 Mpixel -- OK.
* Beat frequency: NOT kHz.  The beat must be moved to >= 4 MHz (most parts) or
  >= 0.75 MHz (epc660).  Then the per-pixel demodulator works, but the 4 phases
  are captured in SEQUENTIAL sub-frames (2 at a time with A/B taps), so any
  motion between sub-frames aliases into phase.
* Photon budget, MLX75027 **[derived]**: 160 ke-/tap x 120 frames/s x 4 raw
  phases ~ 7.7e7 e-/s per pixel -> sigma_phi ~ 1.1e-4 rad/sqrt(s) ~ **6 pm surface
  per sqrt(s)** if shot-limited; the demonstrated system was 8x above shot
  noise, i.e. ~50 pm-class.
* Continuous: yes at 60-120 depth frames/s, MIPI.
* Verdict: **plausible 10-pm-class, not 1-pm-class, and requires relocating the
  beat to MHz** (which costs two AOMs and brings back the sequential-phase
  problem heterodyne was meant to remove).  Key limits: electronics noise, tap
  full well x frame rate, sequential phases.
* Credibility 5-10 y: very high for the sensors existing (Sony, Melexis, ADI,
  Infineon, e2v, ESPROS all ship; automotive and consumer volumes); low for the
  specific needs (no vendor has a reason to add a sub-MHz demodulation mode or
  simultaneous 4-phase capture at 1 Mpixel -- Hydra3D's 3-tap is the nearest).

---

## 3. Coherent / heterodyne lidar focal-plane arrays with per-pixel I/Q or phase

### 3a. Pointcloud GmbH silicon-photonic FMCW FPA (Nature, March 2026)

| Item | Value | Source |
|---|---|---|
| Array | 352 x 176 pixels; >0.6 million photonic components on-chip; 53 um pitch | [Nature 651, 364 (2026)](https://www.nature.com/articles/s41586-026-10183-6); [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Per-pixel detector / electronics | balanced germanium PIN photodetectors in each pixel; pixel-integrated TIA; 50/50 coupler heterodyne mixing; monostatic (Tx and Rx share the pixel) | [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Readout | up to 128 electrical outputs read simultaneously; demonstrated with 8 channels in parallel, rows of 8 pixels read consecutively; ~130 us pixel readout time | [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Frame rate | 3-15 fps possible; 1 fps demonstrated with 8-channel readout | [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Chirp / wavelength / power | 32 us up / 32 us down, 6 GHz total chirp bandwidth; 1310 nm; 32 mW per optical channel, 178 uW on-target per pixel | [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Process | GlobalFoundries 45SPCLO 300 mm silicon photonics with 45 nm SOI RF CMOS | [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Lineage | 512-pixel coherent receiver array (2020, arXiv 2008.02411); 128 x 128 predecessor; "five times increase in pixel count" | [arXiv 2008.02411](https://arxiv.org/pdf/2008.02411); [PMC12979198](https://pmc.ncbi.nlm.nih.gov/articles/PMC12979198/) |
| Status | research demonstration (company: Pointcloud GmbH); not a product | [ISW blog](http://image-sensors-world.blogspot.com/2026/03/pointclouds-4d-fmcw-lidar-paper.html) |

Verdict: this is the most advanced "per-pixel coherent receiver" array in
existence -- but every pixel is a TIME-MULTIPLEXED coherent pixel (<=128 of
62k read at once, 130 us each), at 1310 nm, with Tx on-chip.  Not continuous,
not visible.  Key limit: simultaneous-output count.  5-10 y: the GF 45SPCLO
platform is commercial and the pixel count has gone 512 -> 16k -> 62k in ~6
years, so a 256 x 256 version is credible; a version reading ALL pixels
continuously at a kHz beat is not on any lidar roadmap -- lidar does not need it.

### 3b. MIT Lincoln Laboratory Geiger-mode APD arrays

| Item | Value | Source |
|---|---|---|
| Array formats (framed lidar ROICs) | 32 x 32 (100 um, 350 nm CMOS) -> 128 x 32 -> 256 x 64 -> 256 x 128 (50 um, 180 nm); 256 x 256 demonstrated; "working to scale to megapixel-class" | [Aull et al., IEEE JSTQE 24(2) 2018](https://www.ll.mit.edu/media/4536) |
| Frame rates | early designs 20 kframes/s; newest "approaching 150 kframes/s for the sparsest returns"; 0.5 W typical at 20 kframes/s | [Aull 2018](https://www.ll.mit.edu/media/4536) |
| Passive-imaging counter ROIC | 256 x 256, 25 um, per-pixel 7-bit counter + overflow bit, 180 nm | [Aull 2018](https://www.ll.mit.edu/media/4536) |
| Fielded | MACHETE airborne lidar: dual 256 x 64 InGaAsP/InP arrays at 1064 nm | [AOSTB, DTIC](https://apps.dtic.mil/sti/pdfs/AD1033558.pdf) |
| Coherent detection with the array | "first experimental demonstration of coherent detection using an array of photon-counting detectors ... Geiger-mode APDs with integrated CMOS timing circuitry ... shot-noise-limited detection in the single-return-photon regime" (Ashcom, Kaushik, Heinrichs, COTA 2006); array size n/s in abstract | [Optica COTA 2006 CWB1](https://opg.optica.org/abstract.cfm?URI=COTA-2006-CWB1) |
| Per-pixel count rate | n/s in sources read; Geiger-mode pixels have a hold-off / re-arm dead time per event and 1 event per frame per pixel in framed mode | [Aull 2018](https://www.ll.mit.edu/media/4536) |

Verdict: the heterodyne-with-GmAPD-array idea is real and Lincoln's; but a
Geiger pixel registers at most one photon per frame (<= 150 kframes/s), so the
per-pixel photon budget is <= ~1e5 /s -> sigma_phi >= 3e-3 rad/sqrt(s) ~ 0.2 nm
surface **[derived]** -- three orders short of pm-class.  Built for photon-starved
lidar, the opposite regime from a bright-LO gauge.  Also government lab /
export-controlled hardware (n/c for the specific ROICs).

### 3c. Commercial coherent (FMCW) lidar

* SiLC Eyeonic: "4-channel FMCW LiDAR chip" (Mini); Vista "leverages 8 parallel
  channels" ([Business Wire 2026-05](https://www.businesswire.com/news/home/20260508656346/en/SiLC-Technologies-Launches-Eyeonic-Vista-Ultra-Long-Range-Vision-System-for-Perimeter-Security); [Inside Unmanned Systems](https://insideunmannedsystems.com/silc-technologies-introduces-eyeonic-vision-system-mini-advancing-precision-lidar-technology/)).
* Aeva Atlas / Atlas Ultra / Atlas Orion: FMCW 4D lidar; channel count n/s on
  product pages ([aeva.com/atlas](https://www.aeva.com/atlas/)).
* Insight LiDAR: "chip-scale ... Digital Coherent LiDAR"; channel count n/s
  ([Laser Focus World](https://www.laserfocusworld.com/test-measurement/test-measurement/article/16566872/insight-lidar-launches-fmcw-chip-scale-lidar-for-autonomous-vehicles)).
* Bridger Photonics: gas-mapping lidar, scanning; not an FPA (no FPA product
  found; n/c).

Verdict: all are few-channel SCANNING receivers (4-8 coherent channels), not
focal-plane arrays.  Irrelevant to a 1e4-element continuous read.

---

## 4. Large-channel-count arrays with per-element continuous readout

### 4a. In-pixel digital ROICs (per-pixel ADC / counter)

**Lincoln Laboratory Digital-pixel FPA (DFPA):** 256 x 256, 30 um pitch, each
pixel "includes a full-dynamic-range analog-to-digital converter as well as
local digital signal processing (DSP) support"; the in-pixel ADC is a
pulse-frequency-modulation ADC with a **bidirectional counter** (Cint ~1 fF);
four 256 x 256 designs fabricated in IBM 90 nm CMOS; demonstrated capabilities
include **"On-chip synchronous detection of a pulsed laser in the presence of
strong spatial clutter"**, "Near background shot-noise limited detection at
low-input currents", "~30 mW power dissipation at 100 Hz frame rate", operation
with SWIR through LWIR detectors; 2.5 Gbit/s edge serializers
([Lincoln TechNote DFPA](https://archive.ll.mit.edu/publications/technotes/TechNote_DFPA.pdf);
[DFPA development, ResearchGate](https://www.researchgate.net/publication/253002748_Digital-pixel_Focal_Plane_Array_Development)).

Why it matters: an up/down counter clocked by the beat reference IS a digital
lock-in pixel with unbounded well depth (the counter wraps / is read out) --
exactly the architecture that removes the full-well limit.  Lincoln has shown
synchronous detection on-chip at 256 x 256.  Key limits: access (Lincoln /
government program hardware; commercial licensees n/c), designed for IR
hybrids (a visible-band Si detector tier would be needed), and the counter
clock rate vs a 100 kHz beat (PFM ADC rates n/s).  5-10 y credibility: the
architecture is 15+ years old and fielded; a commercial visible-band DFPA with
a user-clocked bidirectional counter is not known to exist (n/c).

### 4b. Brute-force: high-speed streaming CMOS camera digitizing the beat

| Camera | 128 x 128 frame rate | Full well | Interface | Source |
|---|---|---|---|---|
| Mikrotron EoSens 1.1MCX12 | 23,821 fps at 128 x 128 ROI | n/s | CoaXPress 2.0 | [1stVision](https://www.1stvision.com/cameras/models/Mikrotron/EoSens%201.1MCX12) |
| Mikrotron EoSens 3CXP | 42,450 fps at 128 x 128 | n/s | CoaXPress | [Concurrent EDA](https://www.concurrenteda.com/index.php/gigasens/eosens-3cxp) |
| Mikrotron EoSens 21CXP2 (different sensor) | -- | 120 ke- | 4 x CXP-12 = 50 Gbit/s | [Concurrent EDA](https://www.concurrenteda.com/index.php/gigasens/eosens-21cxp2); [Automation.com](https://www.automation.com/en-us/products/september-2021/mikrotron-21-megapixel-camera-coaxpress-imaging) |

**Fitness [derived]:** sample a 1-5 kHz beat at >= 4 samples/period: a 128 x 128
ROI at 20-40 kfps streams 128 x 128 x 10 bit x 24 kfps = 3.9 Gbit/s -- within
CXP-12, continuous.  Throughput ~ (20-100 ke-) x 2e4 fps = 4e8-2e9 e-/s ->
sigma_phi ~ 2-5e-5 rad/sqrt(s) ~ **1-2.5 pm surface per sqrt(s)** shot-limited.
BUT read noise is paid on every sample (high-speed CMOS read noise n/s here,
typically tens of e-) and the demodulation is done in software from a ~1e5
samples/s/pixel stream.  Pixel count is capped by the fps x pixels product:
512 x 512 at 20 kfps = 5e9 px/s is beyond any CXP camera found.  Verdict:
**exists today for 128 x 128, pm-class on paper, not scalable to 512 x 512;
no engineering risk beyond data handling.**  This is the fallback.

### 4c. Channel-count precedents from other fields

| Precedent | Channels / pixels | Per-channel readout | Source |
|---|---|---|---|
| Butterfly iQ ultrasound-on-chip | 9000-element CMUT array on one CMOS die; "each CMUT is driven by a block with CMOS logic, analog reading circuit and LDMOS transistors"; receiver "analog circuitry supporting each group of eight elements" | continuous RF analog front end (MHz ultrasound) | [Yole flyer](https://medias.yolegroup.com/uploads/2020/04/SP20545-Yole-Butterfly-Network-iQ-CMUT-Ultrasonic-Sensor-Flyer.pdf); [Rothberg et al., PNAS 118(27) 2021](https://www.pnas.org/doi/pdf/10.1073/pnas.2019339118) |
| Pitch-matched 3D-ultrasound ASIC (Delft) | "0.91 mW/Element Pitch-Matched Front-End ASIC with Integrated Subarray Beamforming ADC" | per-element AFE + subarray ADC | [ResearchGate](https://www.researchgate.net/publication/323655693_A_091mWElement_Pitch-Matched_Front-End_ASIC_with_Integrated_Subarray_Beamforming_ADC_for_Miniature_3D_Ultrasound_Probes) |
| CERN Timepix4 (HEP / photon science) | 448 x 512 = 229k pixels, 55 um, 65 nm CMOS; per-pixel ToA LSB 195 ps + ToT; readout up to 163.84 Gbps (16 x 10.24 Gbps); sustains ~3.6 Mhits/mm^2/s; data-driven continuous | per-pixel discriminator + TDC, continuous data-driven | [CERN KT brief](https://knowledgetransfer.web.cern.ch/sites/default/files/technology/timepix4/tech-brief/tech-brieftimepix4.pdf); [Medipix4 page](https://medipix.web.cern.ch/medipix4/); [arXiv 2404.07120](https://arxiv.org/pdf/2404.07120) |

What these prove: 1e4 channels of continuous analog front end on one die
(Butterfly, consumer-priced probe) and 2e5 channels of continuous per-pixel
digital readout at 164 Gbit/s (Timepix4, a collaboration ASIC) are ENGINEERED
and SHIPPED.  Costs: Butterfly probe retail price and Timepix4 per-chip price
n/s / n/c from primary sources; the design-and-mask cost of a 65 nm ASIC of
Timepix4 class is a multi-year, multi-institution effort (CERN + Nikhef + IFAE
per the KT brief).  A custom visible-band "lock-in Timepix" is therefore
credible as a 5-10 y collaboration project, not a purchase.

### 4d. 2-D PIN/APD arrays with per-element TIA+ADC at >= 10 kHz, >= 32 x 32, COTS

None found.  Flash-lidar linear-mode APD ROICs (not surveyed in depth) have
per-pixel bandwidth >> 10 kHz but capture a short range gate per frame, not a
continuous stream; the only continuous per-pixel TIA arrays found are the
silicon-photonic coherent FPA above (time-multiplexed outputs) and the DFPA
(in-pixel ADC).  Flag: absence of evidence, not evidence of absence -- a
targeted vendor query (Hamamatsu, First Sensor, Teledyne, Allegro/Voxtel)
would be needed.

---

## 5. SPAD / photon-counting arrays with per-pixel timing

| Sensor | Pixels / pitch | Timing | Count-rate ceiling per pixel | Frame rates | Status | Source |
|---|---|---|---|---|---|---|
| Pi Imaging SPAD512 | 512 x 512; 16.38 um; PDP 50 % @520 nm; median DCR 25 cps | gate width 6 ns min, 17 ps gate shift | **100 kcps without exposure modulation (240 kcps linearized); 167 Mcps with exposure modulation (400 Mcps linearized)** | 100,000 fps 1-bit for 1 s (buffered); 5,000 fps 4-bit semi-continuous; 400 fps 8-bit continuous | shipping since 2021; USB3 x2 | [SPAD512 datasheet v2.4 (2026)](https://piimaging.com/wp-content/uploads/2024/10/SPAD512Sdatasheet.pdf) |
| Canon 3.2 MP SPAD | 3.2 MP total / 2.1 MP effective; 6.39 um; stacked BSI dToF | 100 ps | n/s | n/s | in MS-500 camera (2023) | [Canon](https://global.canon/en/technology/spad-sensor-2023.html); [TechInsights](https://www.techinsights.com/blog/canon-32-mp-stacked-spad-image-sensor-device-essentials-plus) |
| Sony IMX459 | 597 x 168 (~100k); 10.08 um; PDE 24 % | ~6 ns response per photon | n/s | n/s (dToF lidar, 300 m / 15 cm) | automotive sampling | [Sony release 2021](https://www.sony-semicon.com/en/news/2021/2021090601.html) |

**Fitness [derived]:** pixel counts are fine; the per-pixel photon budget is
not.  SPAD512 continuous: <= 1e5-2.4e5 counts/s/pixel -> sigma_phi >= 2e-3
rad/sqrt(s) ~ **0.1-0.2 nm surface per sqrt(s)**; even at the 167-400 Mcps
exposure-modulated ceiling (a 1 s buffered burst, not continuous) the dead-time
regime gives at best ~1e-4 rad -> ~5 pm, and the readout is 1-bit frames.
Per-pixel time-tagging at 100 ps (Canon) is irrelevant to a kHz beat; what a
heterodyne SPAD gauge would need is an on-chip per-pixel I/Q ACCUMULATOR keyed
to the beat (a "SPAD lock-in"), which exists only as FLIM phasor processing
off-chip.  Verdict: **wrong photon regime** (single-photon detectors with a
bright local oscillator); not a pm path within 5-10 y unless someone builds a
counting lock-in pixel with >1e8 cps linear range -- no such roadmap found.

---

## 6. Ranked short list

| Rank | Path | Fit today | Key limit | 5-10 y credibility | What must be true |
|---|---|---|---|---|---|
| 1 | **Heliotis heliCam C4 / C4M in-pixel analog lock-in** | Only COTS device meeting pixel count (542x512 / 1024x1102), beat band (1-250 kHz), continuous in-pixel integration, external reference input | 10-bit X/Y output (0.1 deg/sample) and n/s read noise; GigE streaming; single vendor; C4 band discrepancy 134 vs 250 kHz | High: 4 generations shipped over >10 years; C4M is a 1 Mpixel lock-in imager | Vendor confirms: (a) long-tc (10-100 ms) operation without output saturation / adjustable gain, (b) measured phase noise per pixel vs photon flux at tc=100 ms is <= 2e-5 rad, (c) 250 kHz applies to the unit sold |
| 2 | **Brute-force streaming high-speed CMOS at 128 x 128 (CXP)** | 20-40 kfps at 128 x 128 exists; 1-2.5 pm/sqrt(s) shot-limited **[derived]** | read noise per sample; software demodulation of ~4 Gbit/s; cannot scale to 512 x 512 | Very high (commodity) | Read noise x sqrt(samples) stays below shot noise; 128 x 128 (1.3 px/actuator) is accepted as "full resolution" |
| 3 | **iToF 2-tap demodulators with the beat moved to 4-100 MHz** (Sony IMX556, Melexis MLX75027, ADI ADSD3100, e2v Hydra3D) | 0.3-1 Mpixel; proven to image optical beats (OPT8241 GW phase camera) | ~6 pm/sqrt(s) shot-limited **[derived]**, electronics-limited 8x worse in the one demonstration; sequential 4-phase capture; needs two AOMs | Very high for the parts; low for a sub-MHz or simultaneous-phase mode | A bright-LO test shows electronics noise << shot noise; sequential-phase aliasing tolerated or 3-tap (Hydra3D) suffices |
| 4 | **Digital-pixel ROIC with in-pixel bidirectional counter (Lincoln DFPA lineage) or a custom "lock-in Timepix"** | 256 x 256 synchronous detection demonstrated by Lincoln; 229k-pixel continuous digital readout demonstrated by Timepix4 | access / ITAR; no visible-band COTS; a custom ASIC is a multi-year collaboration | Medium: architectures proven, product absent | A partner (Lincoln, CERN Medipix, or a ROIC house) takes it on; funding for a 65-90 nm ASIC run |
| 5 | **Silicon-photonic coherent FPA (Pointcloud lineage)** | 352 x 176 pixels with per-pixel balanced PD + TIA, but <=128 simultaneous outputs, 1310 nm | time-multiplexed readout; wavelength; lidar roadmap does not want all-pixel continuous output | Medium for pixel count growth; low for this use | A visible/NIR variant with all-pixel parallel I/Q outputs -- nothing found pointing that way |

Not ranked (unfit): SPAD arrays (count-rate ceiling; wrong photon regime),
Geiger-mode APD coherent arrays (same), commercial FMCW lidar (4-8 channels,
scanning).

---

## 7. Flags -- things I could not confirm

1. **heliCam C4 demodulation band:** datasheet 1-250 kHz vs vendor web page
   305 Hz-134 kHz (and C4M 305 Hz-50 kHz).  Possibly a firmware/model change
   after the 2023 launch.  Ask Heliotis which applies to current stock.
2. **heliCam C4 read noise and the meaning of "FWC 500 ke-"** (per internal
   1-MSa/s sample or per output sample): n/s.  The pm-class conclusion holds
   either way, but the actual noise floor at tc = 100 ms is the number that
   decides rank 1.
3. **heliCam price:** quote only; no public price class found.
4. **ADSD3100 modulation-frequency range:** the ADI datasheet fetch timed out
   twice; left n/c.  ADTF3175 module price $197 @1k (2022) is from a blog
   summary of the ADI release.
5. **ESPROS epc660 MODCLK:** "external LED/LD modulation input" -- whether it
   also drives the pixel demodulation phase is n/c.
6. **Lincoln DFPA commercial availability / visible-band version / counter
   clock rate:** n/c.
7. **Ashcom 2006 coherent GmAPD array size:** n/s in the abstract.
8. **Papers-per-year counts** for Heliotis and for lock-in pixels in general:
   not tallied; only the vendor's paper-list page is cited.
9. **High-speed CMOS read noise and the full well of the specific 128 x 128-ROI
   models:** n/s; the 120 ke- figure is for a different Mikrotron sensor
   (21CXP2).
10. **No independent second vendor of a kHz in-pixel lock-in imager was found**
    -- this may be a search gap rather than a market fact.

---

## Sources (all URLs appear inline above)

Heliotis C4 datasheet; heliotis.com/en/c4; inVISION news; findlight; PMC4957977
(heliCam C3 heterodyne holography); Swissphotonics Heliotis presentation;
arXiv 1907.05224 (optical lock-in camera, GW); arXiv 2105.07886 (OPT8241 phase
camera); Sony IMX556 flyer + iToF page; Melexis MLX75027 datasheet; ADI
ADSD3100 pages; Infineon REAL3 briefs; pmdtec flexx2; Teledyne Hydra3D; ESPROS
epc660 datasheet; Nature s41586-026-10183-6 / PMC12979198 (Pointcloud FPA);
arXiv 2008.02411; Aull et al. IEEE JSTQE 2018 (ll.mit.edu/media/4536); DTIC
AD1033558 (MACHETE); Optica COTA 2006 CWB1; Lincoln TechNote DFPA; Mikrotron
EoSens pages; Yole Butterfly flyer; Rothberg PNAS 2021; CERN Timepix4 KT brief;
PI Imaging SPAD512 datasheet v2.4; Canon SPAD page; Sony IMX459 release; SiLC,
Aeva, Insight LiDAR product/press pages.
