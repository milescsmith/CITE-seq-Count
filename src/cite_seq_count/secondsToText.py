# Gist found here: https://gist.github.com/Highstaker/280a09591df4a5fb1363b0bbaf858f0d


def pluralizeRussian(
    number: float,
    nom_sing: str,
    gen_sing: str,
    gen_pl: str
    ) -> str:
    """
    Changes the hours, minutes, seconds to plural
    """
    s_last_digit = str(number)[-1]

    if int(str(number)[-2:]) in range(11, 20):
        # 11-19
        return gen_pl
    elif s_last_digit == "1":
        # 1
        return nom_sing
    elif int(s_last_digit) in range(2, 5):
        # 2,3,4
        return gen_sing
    else:
        # 5,6,7,8,9,0
        return gen_pl


def secondsToText(secs: float, lang: str = "EN") -> str:
    """
    Converts datetime to human readable hours, minutes, secondes format.

    Args:
            secs (float): Secondes
            lang (string): Language

    Returns:
            string: Human readable datetime format.
    """
    days = secs // 86400
    hours = (secs - days * 86400) // 3600
    minutes = (secs - days * 86400 - hours * 3600) // 60
    seconds = secs - days * 86400 - hours * 3600 - minutes * 60

    if lang == "ES":
        days_text = f'día{"s" if days != 1 else ""}'
        hours_text = f'hora{"s" if hours != 1 else ""}'
        minutes_text = f'minuto{"s" if minutes != 1 else ""}'
        seconds_text = f'segundo{"s" if seconds != 1 else ""}'
    elif lang == "DE":
        days_text = f'Tag{"e" if days != 1 else ""}'
        hours_text = f'Stunde{"n" if hours != 1 else ""}'
        minutes_text = f'Minute{"n" if minutes != 1 else ""}'
        seconds_text = f'Sekunde{"n" if seconds != 1 else ""}'
    elif lang == "RU":
        days_text = pluralizeRussian(days, "день", "дня", "дней")
        hours_text = pluralizeRussian(hours, "час", "часа", "часов")
        minutes_text = pluralizeRussian(minutes, "минута", "минуты", "минут")
        seconds_text = pluralizeRussian(seconds, "секунда", "секунды", "секунд")
    else:
        # Default to English
        days_text = f'day{"s" if days != 1 else ""}'
        hours_text = f'hour{"s" if hours != 1 else ""}'
        minutes_text = f'minute{"s" if minutes != 1 else ""}'
        seconds_text = f'second{"s" if seconds != 1 else ""}'

    return ", ".join(
        filter(
            lambda x: bool(x),
            [
                f"{days} {days_text}" if days else "",
                f"{hours} {hours_text}" if hours else "",
                f"{minutes} {minutes_text}" if minutes else "",
                f"{seconds:.4} {seconds_text}" if seconds else "",
            ],
        )
    )
