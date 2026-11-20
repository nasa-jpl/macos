#!/usr/bin/env bash
# the loop's fiber rows at the CORRECTED budget (3.5e-3 rad/frame) and 3x, baseline camera, both rigs (2026-10-12)
cd "$(dirname "${BASH_SOURCE[0]}")/.."
F="'cophi.loop.cams',{'baseline'},'cophi.loop.rows',{'fiber1','fiber3'},'cophi.loop.fiber',3.5e-3"
./cophi_batch.sh loopfib_lens_plain "'bench.polarizing',false,'stages',{'loop'},$F"
./cophi_batch.sh loopfib_oap_plain  "'bench.optics','oap','bench.polarizing',false,'stages',{'loop'},$F"
echo LOOPFIB_DONE
