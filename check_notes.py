#!/usr/bin/env python3
"""Check that no speaker note is cut off in a build that typesets them.

A beamer note page is one vbox that does not break.  Anything past the
bottom is pushed off the page and never reaches the PDF, and TeX only
warns about it when the overflow is large, so the build log is not a
reliable detector -- three notes were being truncated with only two
warnings to show for it.  This compares the end of every note in the
source against the text of the page it should appear on.

Only the part of the page that holds the note is read.  Several notes
end on the same words as the slide beside them, and reading the whole
page would let the slide vouch for a note that is not there.

Usage: check_notes.py <tex> <pdf> [two-screen|print]

  two-screen  \\NOTES: one page per frame, note in the right half
  print       \\NOTESPRINT: only frames with a note, note in the right
              column of the layout set in the preamble
"""
import re
import subprocess
import sys

TAIL_WORDS = 5
# Left edge of the note, as a fraction of the page width, per layout.
NOTE_LEFT = {'two-screen': 0.50, 'print': 0.48}


def balanced(text, start):
    """Return the contents of the {...} group beginning at text[start]."""
    depth, i = 0, start
    while i < len(text):
        if text[i] == '{' and text[i - 1] != '\\':
            depth += 1
        elif text[i] == '}' and text[i - 1] != '\\':
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
        i += 1
    raise ValueError('unbalanced brace')


def words(s):
    """Letters and digits only, lowercased -- immune to how TeX set it."""
    s = re.sub(r'\\[a-zA-Z]+', ' ', s)
    return re.findall(r'[a-z0-9]+', s.lower())


def notes_by_frame(tex):
    """[(frame number, note text)], frames numbered from 1 as pages are."""
    out, frame = [], 0
    for m in re.finditer(r'\\begin\{frame\}|\\note\{', tex):
        if m.group().startswith(r'\begin'):
            frame += 1
        else:
            out.append((frame, balanced(tex, m.end() - 1)))
    return out


def main(tex_path, pdf_path, layout='two-screen'):
    tex = open(tex_path, encoding='utf-8').read()
    info = subprocess.run(['pdfinfo', pdf_path], capture_output=True,
                          text=True).stdout
    w, h = (float(x) for x in
            re.search(r'Page size:\s+([\d.]+) x ([\d.]+)', info).groups())
    left = int(w * NOTE_LEFT[layout])
    notes = notes_by_frame(tex)
    bad = 0
    for i, (frame, note) in enumerate(notes, start=1):
        # The print build drops the frames that carry no note, so there
        # the n-th note is on page n; the two-screen build keeps them.
        page = i if layout == 'print' else frame
        tail = words(note)[-TAIL_WORDS:]
        text = subprocess.run(
            ['pdftotext', '-q', '-f', str(page), '-l', str(page),
             '-x', str(left), '-y', '0', '-W', str(int(w) - left),
             '-H', str(int(h) + 1), pdf_path, '-'],
            capture_output=True, text=True).stdout
        if ''.join(tail) not in ''.join(words(text)):
            bad += 1
            print(f'  page {page}: note is CUT OFF, missing "{" ".join(tail)}"')
    if bad:
        print(f'{bad} note(s) truncated -- shorten them, the page cannot grow.')
        return 1
    print(f'{layout}: all {len(notes)} notes reach the page intact')
    return 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
