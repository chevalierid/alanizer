# alanizer

Synchronising RTCs to real time:

1. Connect your RTC to your Arduino controller.
2. Download `rtc_set.ino` and `get_time_v3.py` from the Utilities folder.
3. Run `rtc_set.ino` without opening the serial monitor.
4. Run `get_time_v3.py` from the terminal.
5. After waiting a few seconds, open the Arduino serial monitor.
6. Monitor the terminal and serial monitor output until they confirm the RTC has been synced.
7. Run `get_time_v3.py`again to confirm that your RTC is providing the correct time.
