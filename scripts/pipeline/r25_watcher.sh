#!/usr/bin/env bash
# R25 Auto-Launcher (Watcher)
# 逻辑:
#   1. 每 5 分钟检查 R24 是否完成 (final_6_r24.json 存在 或 r24 进程结束)
#   2. R24 完成后, 确认 GPU 空闲, 自动启动 R25
#   3. 防止重复启动 (R25 lock 文件)
#   4. 时间保护: 不早于 06:20 启动 R25 (避免提前抢占, 用户要求)

LOG=/tmp/r25_watcher.log
R24_FINAL=/root/autodl-tmp/r24/final_6_r24.json
R25_LOCK=/tmp/r25_launched.lock
EARLIEST="06:20"

echo "[$(date)] R25 watcher started. Waiting for R24 to finish..." >> $LOG

while true; do
    # 已经启动过 R25 则退出
    if [ -f "$R25_LOCK" ]; then
        echo "[$(date)] R25 already launched. Watcher exiting." >> $LOG
        break
    fi

    # 检查 R24 是否还在跑
    R24_RUNNING=$(ps aux | grep r24_server.py | grep -v grep | wc -l)

    # 检查 R24 是否完成 (final json 存在)
    R24_DONE=0
    if [ -f "$R24_FINAL" ]; then
        R24_DONE=1
    fi
    # R24 进程结束也算完成 (即使没 final, 用 progress)
    if [ "$R24_RUNNING" -eq 0 ]; then
        R24_DONE=1
    fi

    if [ "$R24_DONE" -eq 1 ]; then
        # 时间保护: 当前时间 >= EARLIEST
        NOW=$(date +%H:%M)
        if [[ "$NOW" < "$EARLIEST" ]]; then
            echo "[$(date)] R24 done but waiting until $EARLIEST (now $NOW)" >> $LOG
            sleep 300
            continue
        fi

        # 确认 GPU 空闲 (避免 R24 残留进程)
        GPU_MEM=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
        if [ "$GPU_MEM" -gt 5000 ]; then
            echo "[$(date)] R24 done but GPU still busy ($GPU_MEM MiB). Waiting..." >> $LOG
            sleep 120
            continue
        fi

        # 启动 R25
        echo "[$(date)] R24 done + GPU free + time OK. Launching R25!" >> $LOG
        touch "$R25_LOCK"
        cd /root/autodl-tmp
        nohup stdbuf -oL -eL python3 -u r25_server.py > /tmp/r25.log 2>&1 &
        echo "[$(date)] R25 launched, PID=$!" >> $LOG
        break
    fi

    echo "[$(date)] R24 still running ($R24_RUNNING proc). Check again in 5 min." >> $LOG
    sleep 300
done

echo "[$(date)] Watcher finished." >> $LOG