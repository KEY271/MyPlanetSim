set shell := ["bash", "-euo", "pipefail", "-c"]

root := justfile_directory()
web_dir := root + "/web"
simulator := root + "/build/dev/my_planet_sim"

# 利用可能なタスクを表示する
default:
	@just --list

# C++ と Web の依存関係を初回セットアップする
setup:
	cmake --preset dev
	cmake --build --preset dev
	cd "{{web_dir}}" && npm ci

# 保存済み dataset を読む viewer を起動する
dev ui_port="5173":
	cd "{{web_dir}}" && npm run dev --workspace @myplanetsim/ui -- --host 127.0.0.1 --port "{{ui_port}}"

# 小格子の地形・期間平均 dataset を生成する
dataset:
	#!/usr/bin/env bash
	if [[ ! -x "{{simulator}}" ]]; then
	  cmake --preset dev
	  cmake --build --preset dev
	fi
	"{{simulator}}" --config "{{root}}/configs/phase15_viewer_rest_n4.cfg" --progress-interval-s 0

# C++・Web のブラウザー不要ゲートを実行する
check:
	cmake --build --preset dev
	ctest --preset dev --output-on-failure
	cmake --build build/dev --target format-check
	cd "{{web_dir}}" && npm run lint && npm run typecheck && npm run test && npm run build
	cd "{{web_dir}}" && MPS_VISUAL_DATASET="{{root}}/build/dev/output/phase15-viewer-rest-n4/viewer" npm run test:cli-dataset --workspace @myplanetsim/protocol
