#!/usr/bin/env python3
import json
import unittest

import cereal.messaging as messaging
from openpilot.common.params import Params
from openpilot.selfdrive.manager.process_config import managed_processes

"""
TODO:
- test language settings for all supported languages
"""

class TestNavd(unittest.TestCase):
  def setUp(self):
    self.params = Params()
    self.sm = messaging.SubMaster(['navRoute', 'navInstruction'])

  def tearDown(self):
    managed_processes['navd'].stop()

  def test_simple(self):
    start = {
      "latitude": 32.7427228,
      "longitude": -117.2321177,
    }
    self.params.put("LastGPSPosition", json.dumps(start))

    end = {
      "latitude": 32.7557004,
      "longitude": -117.268002,
    }
    self.params.put("NavDestination", json.dumps(end))

    managed_processes['navd'].start()

    for _ in range(10):
      self.sm.update(1000)
      if all(f > 0 for f in self.sm.rcv_frame.values()):
        break
    else:
      raise Exception("didn't get a route")


if __name__ == "__main__":
  unittest.main()
