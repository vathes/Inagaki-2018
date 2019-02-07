FROM datajoint/jupyter:python3.6

ADD . /src/template_project

RUN pip install -e /src/template_project

