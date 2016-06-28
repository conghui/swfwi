#!/bin/bash
echo "Merging jpg to mp4, waiting..."
Ratio=5
ffmpeg -r ${Ratio} -f image2 -i pic/%03d.jpg -y -b 12M -f mp4 -r ${Ratio} -pix_fmt rgb24 out.mp4
#ffmpeg -r ${Ratio} -f image2 -i pic/proj%03d.jpg -y -s 928x1188 -b 1M -f mp4 -r ${Ratio} -pix_fmt rgb24 out.mp4
