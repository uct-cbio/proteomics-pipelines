#!/usr/bin/env python3

import mqparse
import sys

txt = sys.argv[1]
output = sys.argv[2]

mq_txt = mqparse.mq_txt(txt)

print(txt)
print(output)
