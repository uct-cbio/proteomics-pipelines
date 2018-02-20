#!/usr/bin/env python3

import mqparse
import sys

config = sys.argv[1]
mq_txt = mqparse.mq_txt(config=config)
