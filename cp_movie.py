#!/usr/bin/env python3
"""
Build a daylight-only movie from CoastCam images.

- Recurses station folders (caco-03|caco-04|caco-05) → camera (c1|c2) → year (2023|2024|2025)
  → daily subfolders (e.g., 144_May.24) → image files (*.jpg/*.jpeg)
- Filters by camera and image type and local daylight window
- Filters by local date range
- Sorts by epoch and writes MP4

Edit the CONFIG section near the top to change:
  - station folder (STATION_DIR_NAME)
  - camera (CAMERA)
  - image type (IMAGE_TYPE)
  - date range (START_DATE_STR, END_DATE_STR)
  - twilight mode (TWILIGHT_MODE)
  - output directory (OUTPUT_DIR)
  - FPS, overlay, offsets, etc.

Dependencies (conda-forge):
    conda install -c conda-forge astral opencv tqdm
"""

from pathlib import Path
import re
from datetime import datetime, timezone, date, timedelta
from zoneinfo import ZoneInfo
from typing import List, Dict, Optional, Tuple

import cv2
from tqdm import tqdm

# Astral for dawn/dusk / sunrise/sunset
try:
    from astral import Observer
    from astral.sun import dawn, dusk, sunrise, sunset
except Exception:
    raise RuntimeError(
        "Missing 'astral'. Install in CRS env: conda install -c conda-forge astral"
    )

# ---------------------------------------------------------------------
# CONFIG — edit these values to customize your run
# ---------------------------------------------------------------------

# Root folder that contains station subfolders 'caco-03', 'caco-04', 'caco-05'
MARCONI_ROOT = Path(r"F:\crs\proj\marconi_imagery")

# Station folder to use (choose one of: "caco-03", "caco-04", "caco-05")
STATION_DIR_NAME = "caco-03"      # <-- change to "caco-03" or "caco-05"

# Camera to use ('c1' or 'c2')
CAMERA = "c2"                      # <-- change to "c1" to use camera 1

# Image type to include: one of {"timex", "bright", "dark", "snap", "var"}
IMAGE_TYPE = "timex"               # <-- change as needed

# Year folders to scan under the camera (normally ["2023","2024","2025"])
YEAR_DIRS = ["2023", "2024", "2025", "2026"]

# Date range to include (local dates in America/New_York).
# Use ISO strings "YYYY-MM-DD". Set to None to include all dates.
START_DATE_STR = "2024-08-15"              # e.g., "2024-12-01"
END_DATE_STR   = "2025-01-25"              # e.g., "2025-11-21"

# Local timezone (for dawn/dusk and timestamp overlay)
TZ_NAME = "America/New_York"

# If True, exclude frames before dawn / after dusk (configured below)
FILTER_DAYLIGHT = True

# Twilight mode: choose one of {"civil", "nautical", "astronomical", "sunrise_sunset"}
TWILIGHT_MODE = "civil"

# Location (Outer Cape Cod; Marconi Beach area)
LAT = 41.914
LON = -69.965

# Optional offsets from dawn/dusk (minutes)
DAWN_OFFSET_MIN = 0                # e.g., start 10 min after civil dawn
DUSK_OFFSET_MIN = 0                # e.g., end 10 min before civil dusk

# Video settings
FPS = 20
DRAW_TIMESTAMP = True
OVERLAY_TZ = TZ_NAME
FONT_SCALE = 0.7
FONT_THICKNESS = 2
FONT_COLOR = (255, 255, 255)  # white
SHADOW_COLOR = (0, 0, 0)      # black shadow
BOTTOM_MARGIN = 30
LEFT_MARGIN = 30

# Output directory:
# - Set a full path like Path(r"F:\crs\proj\marconi_imagery\movies") to override.
# - Or set to None to auto-place under station/camera: MARCONI_ROOT/STATION_DIR_NAME/CAMERA/movies
OUTPUT_DIR: Optional[Path] = None

# ---------------------------------------------------------------------
# Regex — matches your filenames exactly, including day-of-month
# Examples:
#   1734714000.Fri.Dec.20_17_00_00.GMT.2024.CACO03.c2.timex.jpg
#   1748131200.Sun.May.25_00_00_00.GMT.2025.CACO04.c1.bright.jpg
# ---------------------------------------------------------------------
NAME_REGEX = re.compile(
    r"""^
    (?P<epoch>\d+)\.                          # epoch at start
    [A-Za-z]{3}\.                             # DOW
    [A-Za-z]{3}\.                             # month
    (?P<day>\d{1,2})_                         # day-of-month
    (?P<hh>\d{2})_(?P<mm>\d{2})_(?P<ss>\d{2}) # HH_MM_SS
    \.GMT\.
    (?P<year>\d{4})\.
    (?P<station>CACO\d{2})\.
    (?P<camera>c\d)\.
    (?P<itype>timex|bright|dark|snap|var)(?:\.[A-Za-z0-9_-]+)?  # allow variant after type
    \.(?:jpe?g)$                              # .jpg/.jpeg (any case via IGNORECASE)
    """,
    re.VERBOSE | re.IGNORECASE
)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def parse_iso_date(dstr: Optional[str]) -> Optional[date]:
    if not dstr:
        return None
    return date.fromisoformat(dstr)

START_DATE = parse_iso_date(START_DATE_STR)
END_DATE   = parse_iso_date(END_DATE_STR)
LOCAL_TZ   = ZoneInfo(TZ_NAME)

def station_code_from_dir(dir_name: str) -> str:
    """
    Convert folder name like 'caco-04' to filename token 'CACO04'
    """
    m = re.match(r"^caco-(\d{2})$", dir_name.lower())
    if not m:
        raise ValueError("STATION_DIR_NAME must be one of 'caco-03', 'caco-04', 'caco-05'")
    return f"CACO{m.group(1)}"

STATION_CODE = station_code_from_dir(STATION_DIR_NAME)

def effective_output_dir() -> Path:
    """
    Return the output directory (either user-specified OUTPUT_DIR or auto path).
    """
    if OUTPUT_DIR is not None:
        if isinstance(OUTPUT_DIR, Path):
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            return OUTPUT_DIR
        else:
            raise TypeError("OUTPUT_DIR must be a Path or None. Use Path(r'...').")
    auto = MARCONI_ROOT / "movies"
    auto.mkdir(parents=True, exist_ok=True)
    return auto

def compose_output_name(station_code: str, camera: str, itype: str,
                        start: Optional[date], end: Optional[date]) -> Path:
    out_dir = effective_output_dir()
    if start and end:
        tag = f"{start.isoformat()}_to_{end.isoformat()}"
    else:
        tag = "all_dates"
    fname = f"{station_code.lower()}_{camera}_{itype}_{tag}_daylight.mp4"
    return out_dir / fname

def find_images_recursive(marconi_root: Path, station_folder: str,
                          camera: str, itype: str,
                          year_dirs: List[str]) -> List[Path]:
    """
    Recursively scan year folders and their subfolders for images that
    contain the station's camera and image type tokens.
    """
    files: List[Path] = []
    camera_token = camera.lower()
    type_token = itype.lower()

    station_root = marconi_root / station_folder
    if not station_root.exists():
        raise FileNotFoundError(f"Station folder does not exist: {station_root}")

    cam_root = station_root / camera_token
    if not cam_root.exists():
        raise FileNotFoundError(f"Camera folder does not exist: {cam_root}")

    for y in year_dirs:
        year_dir = cam_root / str(y)
        if not year_dir.exists():
            print(f"WARNING: Year folder missing: {year_dir}")
            continue

        # Scan recursively for jpg/jpeg (both cases)
        for ext in ("*.jpg", "*.JPG", "*.jpeg", "*.JPEG"):
            for f in year_dir.rglob(ext):
                n = f.name.lower()
                # Quick substring filter to avoid parsing unrelated images
                if camera_token in n and type_token in n:
                    files.append(f)

    # Deduplicate while preserving order
    files = list(dict.fromkeys(files))
    return files

def parse_filename_info(path: Path) -> Optional[Dict]:
    """
    Parse filename to extract epoch, station, camera, type, etc.
    Returns dict or None if not matched.
    """
    m = NAME_REGEX.match(path.name)
    if not m:
        return None
    gd = m.groupdict()

    # Enforce desired camera and image type and station code
    if gd["camera"].lower() != CAMERA.lower():
        return None
    if gd["itype"].lower() != IMAGE_TYPE.lower():
        return None
    if gd["station"].upper() != STATION_CODE.upper():
        return None

    try:
        epoch = int(gd["epoch"])
    except Exception:
        return None

    dt_utc = datetime.fromtimestamp(epoch, tz=timezone.utc)
    return {
        "path": path,
        "epoch": epoch,
        "dt_utc": dt_utc,
        "station": gd["station"].upper(),
        "year": int(gd["year"]),
        "camera": gd["camera"].lower(),
        "itype": gd["itype"].lower(),
    }

def within_date_range(dt_local: datetime,
                      start: Optional[date], end: Optional[date]) -> bool:
    """
    Check if local datetime's date falls within [start, end] inclusive.
    If start or end is None, treat as unbounded.
    """
    d = dt_local.date()
    if start and d < start:
        return False
    if end and d > end:
        return False
    return True

def get_daylight_window(local_date: date, lat: float, lon: float,
                        tzname: str, mode: str) -> Tuple[datetime, datetime]:
    """
    Returns (start_utc, end_utc) for the daylight/twilight window on 'local_date'.
        - "civil": civil dawn to civil dusk (sun at -6°)
        - "nautical": nautical dawn/dusk (-12°)
        - "astronomical": astronomical dawn/dusk (-18°)
        - "sunrise_sunset": sunrise to sunset (sun at 0°)
    Results returned in UTC for comparison with dt_utc in filenames.
    """
    tz = ZoneInfo(tzname)
    obs = Observer(latitude=lat, longitude=lon)

    if mode == "civil":
        d_local = dawn(observer=obs, date=local_date, tzinfo=tz, depression=6)
        s_local = dusk(observer=obs, date=local_date, tzinfo=tz, depression=6)
    elif mode == "nautical":
        d_local = dawn(observer=obs, date=local_date, tzinfo=tz, depression=12)
        s_local = dusk(observer=obs, date=local_date, tzinfo=tz, depression=12)
    elif mode == "astronomical":
        d_local = dawn(observer=obs, date=local_date, tzinfo=tz, depression=18)
        s_local = dusk(observer=obs, date=local_date, tzinfo=tz, depression=18)
    elif mode == "sunrise_sunset":
        d_local = sunrise(observer=obs, date=local_date, tzinfo=tz)
        s_local = sunset(observer=obs, date=local_date, tzinfo=tz)
    else:
        raise ValueError(f"Unknown TWILIGHT_MODE: {mode}")

    # Optional offsets (in local time)
    if DAWN_OFFSET_MIN:
        d_local = d_local + timedelta(minutes=DAWN_OFFSET_MIN)
    if DUSK_OFFSET_MIN:
        s_local = s_local - timedelta(minutes=DUSK_OFFSET_MIN)

    # Convert to UTC
    return (d_local.astimezone(timezone.utc), s_local.astimezone(timezone.utc))

def put_text_with_shadow(img, text, org,
                         font_face=cv2.FONT_HERSHEY_SIMPLEX,
                         font_scale=0.7, color=(255,255,255), thickness=2,
                         shadow_color=(0,0,0)):
    x, y = org
    cv2.putText(img, text, (x+1, y+1), font_face, font_scale,
                shadow_color, thickness+1, cv2.LINE_AA)
    cv2.putText(img, text, (x, y), font_face, font_scale,
                color, thickness, cv2.LINE_AA)

def choose_video_writer(out_path: Path, fps: int, frame_size: Tuple[int, int]) -> cv2.VideoWriter:
    """
    Try H.264/avc1 first (better compression), fall back to mp4v.
    """
    for cc in ('avc1', 'H264', 'mp4v'):
        fourcc = cv2.VideoWriter_fourcc(*cc)
        vw = cv2.VideoWriter(str(out_path), fourcc, fps, frame_size)
        if vw.isOpened():
            print(f"Using codec: {cc}")
            return vw
    raise RuntimeError("Failed to open VideoWriter with avc1/H264/mp4v. "
                       "Ensure OpenCV has FFmpeg/H.264 enabled (conda-forge opencv + ffmpeg).")

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    DRAW_TIMESTAMP=True
    # Validate station/camera folders exist
    station_root = MARCONI_ROOT / STATION_DIR_NAME
    if not station_root.exists():
        raise FileNotFoundError(f"Station folder not found: {station_root}")
    cam_root = station_root / CAMERA
    if not cam_root.exists():
        raise FileNotFoundError(f"Camera folder not found: {cam_root}")

    # 1) Collect candidate files
    files = find_images_recursive(MARCONI_ROOT, STATION_DIR_NAME, CAMERA, IMAGE_TYPE, YEAR_DIRS)
    if not files:
        print(f"No {IMAGE_TYPE} images for {STATION_DIR_NAME}/{CAMERA} found.")
        return

    print(f"Found {len(files)} candidate {STATION_DIR_NAME}/{CAMERA}/{IMAGE_TYPE} images (recursive).")

    # 2) Parse valid files via regex and station/camera/type enforcement
    parsed: List[Dict] = []
    failed = 0
    for f in files:
        info = parse_filename_info(f)
        if info:
            parsed.append(info)
        else:
            failed += 1

    if not parsed:
        print(f"No files matched the expected naming pattern for {STATION_DIR_NAME}/{CAMERA}/{IMAGE_TYPE}.")
        print(f"Regex failures: {failed}")
        return

    print(f"Parsed {len(parsed)} files; {failed} skipped due to non-matching names.")

    # 3) Sort by epoch (ascending)
    parsed.sort(key=lambda x: x["epoch"])

    # 4) Filter by date range and daylight (optional)
    tz = LOCAL_TZ
    daylight_cache: Dict[date, Tuple[datetime, datetime]] = {}

    included: List[Dict] = []
    excluded_night = 0
    excluded_date = 0

    for item in parsed:
        dt_local = item["dt_utc"].astimezone(tz)
        if not within_date_range(dt_local, START_DATE, END_DATE):
            excluded_date += 1
            continue

        if FILTER_DAYLIGHT:
            local_d = dt_local.date()
            if local_d not in daylight_cache:
                try:
                    start_utc, end_utc = get_daylight_window(local_d, LAT, LON, TZ_NAME, TWILIGHT_MODE)
                except Exception:
                    start_utc, end_utc = None, None
                daylight_cache[local_d] = (start_utc, end_utc)

            start_utc, end_utc = daylight_cache[local_d]
            if start_utc is None or end_utc is None:
                excluded_night += 1
                continue

            if not (start_utc <= item["dt_utc"] <= end_utc):
                excluded_night += 1
                continue

        included.append(item)

    if not included:
        print("No images within the requested date range and daylight window.")
        print(f"Excluded (date): {excluded_date}, Excluded (night): {excluded_night}")
        return

    print(f"Included frames: {len(included)} "
          f"(excluded date: {excluded_date}, excluded night: {excluded_night})")

    # 5) Initialize video writer
    first_img = cv2.imread(str(included[0]["path"]))
    if first_img is None:
        print(f"Failed to read first image: {included[0]['path']}")
        return
    height, width = first_img.shape[:2]

    out_path = compose_output_name(STATION_CODE, CAMERA, IMAGE_TYPE, START_DATE, END_DATE)
    vw = choose_video_writer(out_path, FPS, (width, height))

    # 6) Write frames
    for item in tqdm(included, desc="Writing video"):
        img = cv2.imread(str(item["path"]))
        if img is None:
            continue

        # Resize if needed
        if img.shape[1] != width or img.shape[0] != height:
            img = cv2.resize(img, (width, height), interpolation=cv2.INTER_AREA)

        if DRAW_TIMESTAMP:
            dt_local = item["dt_utc"].astimezone(ZoneInfo(OVERLAY_TZ))
            ts = dt_local.strftime("%Y-%m-%d %H:%M:%S %Z")
            put_text_with_shadow(
                img, ts,
                org=(LEFT_MARGIN, height - BOTTOM_MARGIN),
                font_scale=FONT_SCALE,
                color=FONT_COLOR,
                thickness=FONT_THICKNESS,
                shadow_color=SHADOW_COLOR
            )

        vw.write(img)

    vw.release()
    print(f"Done. Wrote: {out_path}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.")