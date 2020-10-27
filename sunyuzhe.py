# author 孙雨哲, 2020.10.27

from sys import argv
from pathlib import Path

print('Usage: python sunyuzhe.py file1 file2 10')
print('File1: forward fastq')
print('File2: reverse fastq')
print('10: barcode length')
disclaimer = (b'\xe6\x9c\xac\xe4\xba\xba\xe6\x89\xbf\xe8\xaf\xba\xef\xbc\x8c'
              b'\xe4\xbd\xbf\xe7\x94\xa8\xe8\xaf\xa5\xe7\xa8\x8b\xe5\xba\x8f'
              b'\xe4\xbf\xae\xe6\x94\xb9\xe6\x95\xb0\xe6\x8d\xae\xe9\x80\xa0'
              b'\xe6\x88\x90\xe7\x9a\x84\xe4\xb8\x80\xe5\x88\x87\xe5\x90\x8e'
              b'\xe6\x9e\x9c\xef\xbc\x8c\xe7\x94\xb1\xe6\x9c\xac\xe4\xba\xba'
              b'\xe6\x89\xbf\xe6\x8b\x85\xef\xbc\x8c\xe4\xb8\x8e\xe7\xa8\x8b'
              b'\xe5\xba\x8f\xe4\xbd\x9c\xe8\x80\x85\xe6\x97\xa0\xe5\x85\xb3'
              b'(\xe8\xbe\x93\xe5\x85\xa5 y \xe7\xbb\xa7\xe7\xbb\xad\xe8\xbf'
              b'\x90\xe8\xa1\x8c\xef\xbc\x89\xef\xbc\x9a')
ok = input(disclaimer.decode('utf8'))
if not ok.upper().startswith('Y'):
    raise SystemExit(0)
barcode_len = int(argv[3])
f = open(argv[1], 'rb')
r = open(argv[2], 'rb')
out_file = Path(argv[2]).with_suffix('.edit.fastq')
edited_r = open(out_file, 'wb')
iter_ = zip(f, r)
for line_f, line_r in iter_:
    # head
    edited_r.write(line_r)
    line_f, line_r = next(iter_)
    edited_line = b''.join([line_f[:barcode_len], line_r[barcode_len:]])
    # seq
    edited_r.write(edited_line)
    line_f, line_r = next(iter_)
    # +
    edited_r.write(line_r)
    line_f, line_r = next(iter_)
    # qual
    edited_line = b''.join([line_f[:barcode_len], line_r[barcode_len:]])
    edited_r.write(line_r)
print('Done.')
print('Output file:', out_file)
