#!/usr/bin/env python3
import os
import re
import subprocess
import sys
from typing import Iterable, List, Optional

from tqdm import tqdm

from openpilot.selfdrive.car.tests.routes import routes as test_car_models_routes
from openpilot.selfdrive.test.helpers import sanitize, apply_to_log
from openpilot.selfdrive.test.process_replay.test_processes import source_segments as replay_segments
from openpilot.tools.lib.azure_container import AzureContainer
from openpilot.tools.lib.openpilotcontainers import DataCIContainer, DataProdContainer, OpenpilotCIContainer, OpenpilotSanetizedDataset

SOURCES: List[AzureContainer] = [
  DataProdContainer,
  DataCIContainer,
  OpenpilotCIContainer
]

DEST = OpenpilotSanetizedDataset

def upload_route(path: str, exclude_patterns: Optional[Iterable[str]] = None) -> None:
  if exclude_patterns is None:
    exclude_patterns = [r'dcamera\.hevc']

  r, n = path.rsplit("--", 1)
  r = '/'.join(r.split('/')[-2:])  # strip out anything extra in the path
  destpath = f"{r}/{n}"
  for file in os.listdir(path):
    if any(re.search(pattern, file) for pattern in exclude_patterns):
      continue
    DEST.upload_file(os.path.join(path, file), f"{destpath}/{file}")


def sync_to_ci_public(route: str) -> bool:
  dest_container, dest_key = DEST.get_client_and_key()
  key_prefix = route.replace('|', '/')
  dongle_id = key_prefix.split('/')[0]

  if next(dest_container.list_blob_names(name_starts_with=key_prefix), None) is not None:
    return True

  print(f"Uploading {route}")
  for source_container in SOURCES:
    # Get all blobs (rlogs) for this route, strip personally identifiable data, and upload to CI
    print(f"Downloading {route}")
    source_key = None
    for source in SOURCES:
      source_container, key = source.get_client_and_key()
      print(f"Trying {source_container.url}")
      blobs = list(source_container.list_blob_names(name_starts_with=key_prefix))
      blobs = [b for b in blobs if re.match(r".*/rlog.bz2", b)]
      print(f"Found {len(blobs)} segments")
      if len(blobs):
        break
    else:
      print("No segments found in source containers")
      print("Failed")
      return False

    for blob_name in blobs:
      if re.search(r"rlog|qlog", blob_name):
        print('downloading', blob_name)
        data = source_container.download_blob(blob_name).readall()
        data = apply_to_log(data, sanitize)

        print(f"Uploading {blob_name} to {dest_container.url}")
        dest_container.upload_blob(blob_name, data)
      else:
        print('copying', blob_name)
        dest_blob_client = dest_container.get_blob_client(blob_name)
        print(source_container.get_blob_client(blob_name).url)
        dest_blob_client.start_copy_from_url(f"{source_container.get_blob_client(blob_name).url}?{source_key}")

    print("Success")
    return True

  return False


if __name__ == "__main__":
  failed_routes = []

  to_sync = sys.argv[1:]

  if not len(to_sync):
    # sync routes from the car tests routes and process replay
    to_sync.extend([rt.route for rt in test_car_models_routes])
    to_sync.extend([s[1].rsplit('--', 1)[0] for s in replay_segments])
  
  for r in tqdm(to_sync):
    if not sync_to_ci_public(r):
      failed_routes.append(r)

  if len(failed_routes):
    print("failed routes:", failed_routes)
