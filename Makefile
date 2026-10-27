PYTHON := .venv/bin/python

.PHONY: setup pipeline dashboard

setup:
	python3 -m venv .venv
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) frequencies.py
	$(PYTHON) responder_analysis.py
	$(PYTHON) subset_analysis.py

dashboard:
	$(PYTHON) -m streamlit run app.py --server.headless true --server.port 8501 --server.address 0.0.0.0
