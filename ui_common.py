import base64

import streamlit as st
from streamlit.source_util import get_pages


def hide_666_pages():
    main_script_path_str = "."
    current_pages = get_pages(main_script_path_str)
    for idx, it in enumerate(current_pages.items()):
        key, value = it
        if "666" in value['script_path']:
            # A good ol' hack
            st.write(
                f"""<style>
                                    li:nth-child({idx + 1}) {{
                                        display: none;
                                        visibility: hidden;
                                    }}
                                </style>
                                """,
                unsafe_allow_html=True,
            )


def on_all_pages():
    hide_666_pages()


def link_svg(link, svg):
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<a href="%s"><img src="data:image/svg+xml;base64,%s"/></a>' % (link, b64)
    st.write(html, unsafe_allow_html=True)
