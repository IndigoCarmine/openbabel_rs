// Language picker for the openbabel-rs guide.
//
// Injects a small dropdown into the mdBook menu bar that switches between the
// English book (served at the site root) and the Japanese book (served under
// `ja/`). It computes the equivalent page in the other language from the current
// URL, so it works both locally (`mdbook serve`) and on GitHub Pages (where the
// site lives under a project path such as `/openbabel_rs/`).
//
// Languages are built to:  <site-base>/<page>        (en, the source language)
//                          <site-base>/ja/<page>     (ja, from po/ja.po)
(function () {
    "use strict";

    var LANGS = [
        { code: "en", label: "English", prefix: "" },
        { code: "ja", label: "日本語", prefix: "ja/" },
    ];

    // Work out, for the current page: the site base URL, this page's path within
    // its book, and which language we are currently viewing.
    function currentContext() {
        var here = new URL(window.location.href);
        // `path_to_root` is a global mdBook sets on every page: the relative path
        // from this page back to its book root.
        var rootRel = (typeof path_to_root !== "undefined" && path_to_root) || "./";
        var bookRoot = new URL(rootRel, here);

        var rel = here.pathname.slice(bookRoot.pathname.length);
        if (rel === "") {
            rel = "index.html";
        }

        var current = "en";
        var base = bookRoot.pathname; // for `en`, the book root IS the site base
        for (var i = 0; i < LANGS.length; i++) {
            var prefix = LANGS[i].prefix;
            if (prefix && bookRoot.pathname.endsWith("/" + prefix)) {
                current = LANGS[i].code;
                base = bookRoot.pathname.slice(0, bookRoot.pathname.length - prefix.length);
                break;
            }
        }
        return { origin: here.origin, base: base, rel: rel, hash: here.hash, current: current };
    }

    function labelFor(code) {
        for (var i = 0; i < LANGS.length; i++) {
            if (LANGS[i].code === code) {
                return LANGS[i].label;
            }
        }
        return code;
    }

    function build() {
        var menuRight = document.querySelector(".menu-bar .right-buttons");
        if (!menuRight) {
            return;
        }
        var ctx = currentContext();

        var wrapper = document.createElement("div");
        wrapper.className = "lang-picker-wrapper";

        var button = document.createElement("button");
        button.className = "icon-button lang-button";
        button.id = "language-toggle";
        button.type = "button";
        button.title = "Change language / 言語を切り替え";
        button.setAttribute("aria-label", "Change language");
        button.setAttribute("aria-haspopup", "menu");
        button.setAttribute("aria-expanded", "false");
        button.textContent = "🌐 " + labelFor(ctx.current);

        var list = document.createElement("ul");
        list.id = "language-list";
        list.className = "theme-popup";
        list.setAttribute("role", "menu");
        list.setAttribute("aria-label", "Languages");

        LANGS.forEach(function (l) {
            var li = document.createElement("li");
            li.className = "theme";
            li.setAttribute("role", "none");

            var a = document.createElement("a");
            a.setAttribute("role", "menuitem");
            a.textContent = l.label;
            a.href = ctx.origin + ctx.base + l.prefix + ctx.rel + ctx.hash;
            if (l.code === ctx.current) {
                li.setAttribute("aria-current", "true");
            }

            li.appendChild(a);
            list.appendChild(li);
        });

        wrapper.appendChild(button);
        wrapper.appendChild(list);
        menuRight.insertBefore(wrapper, menuRight.firstChild);

        function close() {
            list.style.display = "none";
            button.setAttribute("aria-expanded", "false");
        }
        function open() {
            list.style.display = "block";
            button.setAttribute("aria-expanded", "true");
        }
        close();

        button.addEventListener("click", function (e) {
            e.preventDefault();
            e.stopPropagation();
            if (list.style.display === "block") {
                close();
            } else {
                open();
            }
        });
        document.addEventListener("click", function (e) {
            if (!wrapper.contains(e.target)) {
                close();
            }
        });
        document.addEventListener("keydown", function (e) {
            if (e.key === "Escape") {
                close();
            }
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", build);
    } else {
        build();
    }
})();
