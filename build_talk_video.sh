#!/bin/sh
# Renders the two films used by the talk and encodes them for beamer.
#
#   talk_8cell_spin.mp4     slide 2, a turntable of the finished net
#   talk_8cell_unfold.mp4   slide 3, the polytope opening into the net
#                           one BFS layer at a time, plus
#                           talk_8cell_unfold_poster.png, the frame that
#                           stands in for it in the plain PDF
#   talk_mainresult_spin.mp4  slide 6, Figure 1 with its two large
#                           panels turning, plus talk_mainresult_poster.png
#
# The frames are throwaway; only the mp4s are kept, and they must sit
# next to the slides, because pympress resolves a \movie filename
# relative to the PDF.
#
# H.264 in MP4 with yuv420p: that is what the GStreamer stack behind
# pympress decodes here (qtdemux + avdec_h264, with VideoToolbox also
# available).  Keep the dimensions even or yuv420p will not encode.
set -e
cd "$(dirname "$0")"

WORK="${TMPDIR:-/tmp}/talkvid_$$"
trap 'rm -rf "$WORK"' EXIT

OUT=600          # final side, in pixels
MARGIN=104       # per cent of the tight crop to keep around the drawing

# An explicit PlotRange in Graphics3D clips, so it cannot be tightened
# past the geometry; and the projection of a box is a hexagon, whose
# bounding box is much larger than the drawing inside it.  A third of
# the picture is therefore empty whatever the plot range.  We recover it
# here: cropdetect over every frame gives the union of what is actually
# drawn, and we crop to the square around that before scaling down.
# Measuring rather than hard-coding means the numbers cannot go stale
# when the camera or the timeline changes.
tight_crop() {   # tight_crop <frame dir> [margin %] -> "side:side:x:y"
  margin=${2:-$MARGIN}
  det=$(ffmpeg -loglevel info -i "$1/f%04d.png" \
          -vf cropdetect=limit=45:round=2:reset=0 -f null - 2>&1 |
        grep -o 'crop=[0-9]*:[0-9]*:[0-9]*:[0-9]*' | tail -1)
  full=$(ffprobe -v error -select_streams v -show_entries stream=width \
           -of csv=p=0 "$1/f0000.png")
  echo "${det#crop=} $full $margin" | awk '{
    split($1, c, ":"); w=c[1]; h=c[2]; x=c[3]; y=c[4]; s=$2;
    cx = x + w/2; cy = y + h/2;
    half = int((w > h ? w : h) / 2 * $3 / 100);
    if (2*half > s) half = int(s/2);
    x0 = int(cx - half); y0 = int(cy - half);
    if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
    if (x0 + 2*half > s) x0 = s - 2*half;
    if (y0 + 2*half > s) y0 = s - 2*half;
    side = 2*half - (2*half) % 2;
    printf "%d:%d:%d:%d", side, side, x0, y0;
  }'
}

encode() {   # encode <frame dir> <output> [crop spec]
  if [ -n "$3" ]; then vf="crop=$3,scale=$OUT:$OUT"; else vf="scale=$OUT:$OUT"; fi
  ffmpeg -y -loglevel error -framerate 30 -i "$1/f%04d.png" \
         -vf "$vf" -c:v libx264 -preset slow -crf 20 -pix_fmt yuv420p \
         -movflags +faststart "$2"
}

echo "--- slide 2: turntable"
wolframscript -file src/figures/make_talk_video_8cell.wls "$WORK/spin"
# Not cropped: this one is framed to match face_rotation_net_8cell.png,
# which is its poster, and cropping would break the match.
encode "$WORK/spin" talk_8cell_spin.mp4

echo "--- slide 3: unfolding"
wolframscript -file src/figures/make_talk_video_unfold.wls "$WORK/unfold"
CROP=$(tight_crop "$WORK/unfold")
echo "    cropping to $CROP before scaling to $OUT"
encode "$WORK/unfold" talk_8cell_unfold.mp4 "$CROP"
# Frame 100 falls in the pause between the two layers: six cells are out
# and flat, the seventh is still out of the hyperplane.  As a still that
# says more than either end of the film, and it is the picture the slide
# shows when the movie cannot be played.
ffmpeg -y -loglevel error -start_number 100 -i "$WORK/unfold/f%04d.png" \
       -frames:v 1 -vf "crop=$CROP,scale=$OUT:$OUT" talk_8cell_unfold_poster.png

echo "--- slide 6: the two large panels of Figure 1"
# The figure is an Illustrator assembly, so the only way to animate part
# of it is to lay new panels back into it.  These rectangles are the two
# lower panels of 260612AllFaceRotation.pdf, measured at 100 dpi on the
# rendered page (1157x874, panels 562x560 at x=12 and x=583, y=302).
# Re-measure them if the figure is ever redrawn.
PW=562; PH=560; PY=302; PX1=12; PX2=583
pdftoppm -r 100 -png -singlefile 260612AllFaceRotation.pdf "$WORK/figure"
wolframscript -file src/figures/make_talk_video_mainresult.wls "$WORK/main" 240
C120=$(tight_crop "$WORK/main/c120" 104)
C600=$(tight_crop "$WORK/main/c600" 104)
echo "    120-cell panel cropped to $C120, 600-cell to $C600"
ffmpeg -y -loglevel error \
  -loop 1 -framerate 30 -i "$WORK/figure.png" \
  -framerate 30 -i "$WORK/main/c120/f%04d.png" \
  -framerate 30 -i "$WORK/main/c600/f%04d.png" \
  -filter_complex "[1]crop=$C120,scale=$PW:$PH[a];\
                   [2]crop=$C600,scale=$PW:$PH[b];\
                   [0][a]overlay=$PX1:$PY:shortest=1[t];\
                   [t][b]overlay=$PX2:$PY,scale=900:680[v]" \
  -map "[v]" -c:v libx264 -preset slow -crf 20 -pix_fmt yuv420p \
  -movflags +faststart talk_mainresult_spin.mp4
# Frame 0 stands where the printed figure stands, so the slide looks
# unchanged until the film is started.
ffmpeg -y -loglevel error -i talk_mainresult_spin.mp4 -frames:v 1 \
       talk_mainresult_poster.png

echo
ls -lh talk_8cell_spin.mp4 talk_8cell_unfold.mp4 talk_8cell_unfold_poster.png \
       talk_mainresult_spin.mp4 talk_mainresult_poster.png
for f in talk_8cell_spin.mp4 talk_8cell_unfold.mp4 talk_mainresult_spin.mp4; do
  echo "$f: $(ffprobe -v error -show_entries stream=width,height,nb_frames \
        -of csv=p=0:s=x "$f")"
done
