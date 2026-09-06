#!/bin/sh
# Builds the three PDFs wanted for the JCDCG^3 talk.
#
#   jcdcg3_slides.pdf        plain 4:3, no notes -- hand this over if the
#                            room supplies the machine
#   jcdcg3_slides_notes.pdf  wide, slide left and note right -- open this
#                            one with pympress on the laptop
#   jcdcg3_slides_print.pdf  each script beside its slide, one page per
#                            script -- print it and take it up with you,
#                            because handing over the plain PDF means
#                            there is no presenter screen to read
#
# Run twice through pdflatex each: metropolis needs the second pass for
# the progress bar and the frame numbering.
set -e
cd "$(dirname "$0")"

run() {   # run <jobname> <preamble-macros>
  for pass in 1 2; do
    # Quiet unless it breaks; then show the error rather than the exit
    # status, since the usual cause is a typo in a \note.
    if ! pdflatex -interaction=nonstopmode -halt-on-error \
                  -jobname="$1" "$2\\input{jcdcg3_slides.tex}" >/dev/null
    then
      echo "LaTeX stopped while building $1.pdf:"
      grep -A4 -m2 '^!' "$1.log"
      echo "(nothing was overwritten; fix the source and run this again)"
      exit 1
    fi
  done
  echo "$1.pdf: $(pdfinfo "$1.pdf" | awk '/^Pages/{print $2" pages"} /^Page size/{print $3"x"$5" pt"}' | paste -sd' ' -)"
}

run jcdcg3_slides       ""
run jcdcg3_slides_notes "\\def\\NOTES{}"
run jcdcg3_slides_print "\\def\\NOTESPRINT{}"

echo
echo "warnings (overfull/underfull, expected: 1 at the title frame):"
grep -c "Overfull\|Underfull" jcdcg3_slides.log || true

echo
# A note page is one vbox and does not break, so a note that is too long
# is pushed off the page and never reaches the PDF.  TeX warns only when
# the overflow is large, so the log misses the near misses; this reads
# the PDF back instead.
python3 check_notes.py jcdcg3_slides.tex jcdcg3_slides_notes.pdf two-screen
python3 check_notes.py jcdcg3_slides.tex jcdcg3_slides_print.pdf print

echo
echo "carry ALL THREE mp4s NEXT TO the pdf --"
echo "the movies on slides 2, 3 and 6 are referenced by relative name, and"
echo "without them the posters just sit there.  ./build_talk_video.sh remakes"
echo "them, together with the two poster PNGs the slides embed."
echo
echo "print jcdcg3_slides_print.pdf (10 pages, landscape, fit to page) and"
echo "carry it: on the room's machine there is no presenter screen."
echo
echo "to present:  pympress -t 12:30 -N right jcdcg3_slides_notes.pdf"
echo "  n cycle notes mode   s swap the two screens   b blank the projector"
echo "  p pause timer   r reset timer   f fullscreen the slide window"
