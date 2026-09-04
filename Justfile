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

# offline/mock UI だけを起動する
ui ui_port="5173":
	cd "{{web_dir}}" && npm run dev --workspace @myplanetsim/ui -- --host 127.0.0.1 --port "{{ui_port}}"

# native gateway だけを起動する（token は標準出力の JSON に表示される）
gateway gateway_port="8787":
	#!/usr/bin/env bash
	if [[ ! -x "{{simulator}}" ]]; then
	  cmake --preset dev
	  cmake --build --preset dev
	fi
	if [[ ! -d "{{web_dir}}/node_modules" ]]; then
	  cd "{{web_dir}}"
	  npm ci
	fi
	cd "{{web_dir}}"
	MPS_SIMULATOR_BINARY="{{simulator}}" \
	MPS_REST_PRESET="{{root}}/configs/phase3_rest_n4.cfg" \
	MPS_DRY_PRESET="{{root}}/configs/phase5_visualizer_rest_n4.cfg" \
	MPS_RUN_ROOT="{{root}}/.runs" \
	MPS_GATEWAY_PORT="{{gateway_port}}" \
	npm run start --workspace @myplanetsim/gateway

# native gateway と live UI を同じ session token で起動する
[no-exit-message]
dev gateway_port="8787" ui_port="5173":
	#!/usr/bin/env bash
	if [[ ! -x "{{simulator}}" ]]; then
	  cmake --preset dev
	  cmake --build --preset dev
	fi
	if [[ ! -d "{{web_dir}}/node_modules" ]]; then
	  cd "{{web_dir}}"
	  npm ci
	fi
	cd "{{web_dir}}"
	token="$(node --input-type=module -e 'import { randomBytes } from "node:crypto"; process.stdout.write(randomBytes(32).toString("hex"))')"
	export MPS_SIMULATOR_BINARY="{{simulator}}"
	export MPS_REST_PRESET="{{root}}/configs/phase3_rest_n4.cfg"
	export MPS_DRY_PRESET="{{root}}/configs/phase5_visualizer_rest_n4.cfg"
	export MPS_RUN_ROOT="{{root}}/.runs"
	export MPS_GATEWAY_PORT="{{gateway_port}}"
	export MPS_SESSION_TOKEN="$token"
	npm run start --workspace @myplanetsim/gateway &
	gateway_pid=$!
	cleanup() {
	  kill "$gateway_pid" 2>/dev/null || true
	  wait "$gateway_pid" 2>/dev/null || true
	}
	on_signal() { exit 0; }
	trap cleanup EXIT
	trap on_signal INT TERM
	ready=false
	for _ in $(seq 1 50); do
	  if ! kill -0 "$gateway_pid" 2>/dev/null; then
	    wait "$gateway_pid"
	  fi
	  if TOKEN="$token" PORT="{{gateway_port}}" node --input-type=module -e 'const response = await fetch(`http://127.0.0.1:${process.env.PORT}/api/v1/capabilities`, { headers: { authorization: `Bearer ${process.env.TOKEN}` } }); process.exit(response.ok ? 0 : 1)' 2>/dev/null; then
	    ready=true
	    break
	  fi
	  sleep 0.1
	done
	if [[ "$ready" != true ]]; then
	  echo "gateway did not become ready" >&2
	  exit 1
	fi
	ui_path="/#token=$token&gateway=http%3A%2F%2F127.0.0.1%3A{{gateway_port}}"
	echo
	echo "Live UI: http://127.0.0.1:{{ui_port}}$ui_path"
	echo
	echo "NOTE: the plain http://127.0.0.1:{{ui_port}}/ that Vite prints below has no session"
	echo "      token and falls back to the offline demo, which serves only the shallow-water"
	echo "      preset. Use the Live UI link above; --open below opens it for you."
	echo "Press Ctrl-C to stop the UI and gateway."
	echo
	npm run dev --workspace @myplanetsim/ui -- --host 127.0.0.1 --port "{{ui_port}}" --open "$ui_path"

# C++・Web・native live の主要ゲートを実行する
check:
	cmake --build --preset dev
	ctest --preset dev --output-on-failure
	cmake --build build/dev --target format-check
	cd "{{web_dir}}" && npm run lint && npm run test && npm run build
	cd "{{web_dir}}" && npm run test:live --workspace @myplanetsim/gateway
	cd "{{web_dir}}" && npm run test:browser --workspace @myplanetsim/ui
